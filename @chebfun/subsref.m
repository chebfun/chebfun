function varargout = subsref(f, index)
%SUBSREF   CHEBFUN subsref.
% ( )
%   F(X) returns the values of the CHEBFUN F evaluated on the array X. If X
%   falls on a breakpoint of F, the corresponding value from F.IMPULSES is
%   returned. F(X, 'left') or F(X, 'right') will evaluate F at breakpoints
%   using left- or right-hand limits, respectively. See CHEBFUN/FEVAL for
%   further details. F(:) returns F.
%
%   If F is an array-valued CHEBFUN then F(X, COL) returns the values of the
%   columns specified by the vector COL at the points X. Similarly, F(:, COL)
%   returns a new array-vlaued CHEBFUN containing only the columns specified in
%   COL. In both cases, COL should be a row vector.
%
%   F(G), where G is also a CHEBFUN, computes the composition of F and G. See
%   CHEBFUN/COMPOSE for further details.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2} restricts F to the domain [S1, S2] < [F.ENDS(1), F.ENDS(end)] and
%   simplifies the result. See RESTRICT  and SIMPLIFY for further details. Note
%   that F{[S1, S2]} is not supported due to the behaviour of the MATLAB
%   subsref() command.
%
% See also FEVAL, COMPOSE, GET, RESTRICT, SIMPLIFY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document for array-valued CHEBFUN objects and quasimatrices.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'
        
        % Deal with row CHEBFUN objects:
        isTransposed = f(1).isTransposed;
        if ( isTransposed )
            f = f.';
            if ( length(idx) > 1 )
                idx(1:2) = idx([2,1]);
            end
        end
        
        % Where to evaluate:
        x = idx{1}; 
        
        % Initialise:
        numCols = numColumns(f);
        columnIndex = 1:numCols; % Column index for array-valued CHEBFUNs.
        varIn = {};              % Additional arguments.
        
        % Deal with additional arguments:
        if ( (length(idx) == 2) && ...
             any(strcmpi(idx{2}, {'left', 'right', '-', '+'})) )
            % f(x, 'left') or f(x, 'right'):
            varIn = {idx(2)};
            
        elseif ( length(idx) == 2 )
            % f(x, m), for array-valued CHEBFUN objects:
            if ( strcmp(idx{2}, ':') )
                % Do nothing, as this has already been done above:
                % columnIndex = 1:size(f, 2);
            elseif ( max(idx{2}) <= columnIndex(end) )
                columnIndex = idx{2};
                if ( size(columnIndex, 2) == 1 )
                    columnIndex = columnIndex.';
                end
                if ( size(columnIndex, 1) > 1 )
                    error('CHEBFUN:CHEBFUN:subsref:colidx', ...
                        'Column index must be a vector of integers.')
                end               
            elseif ( max(idx{2}) > columnIndex(end) )
                error('CHEBFUN:CHEBFUN:subsref:badsubscript', ...
                    'Index exceeds CHEBFUN dimensions.');
            elseif ( isempty(idx{2}) )
                columnIndex = [];
            end

        elseif ( length(idx) == 2 && strcmp(idx{2}, ':') )
            % This is OK.
            
        elseif ( length(idx) > 1 )
            error('CHEBFUN:CHEBFUN:subsref:dimensions', ...
                'Index exceeds CHEBFUN dimensions.')
            
        end

        % Compute the output:
        if ( isnumeric(x) )
            % Call FEVAL():
            out = feval(f, x, varIn{:});

            % Figure out which columns of the output we need to select:
            % (NB:  This code uses the assumption that columnIndex is a row.)
            evalPointCols = size(x, 2);
            outCols = bsxfun(@plus, (columnIndex - 1)*evalPointCols + 1, ...
                (0:1:(evalPointCols - 1)).');

            % Select only the required columns of the output:
            % (NB:  The cell array of colons is a hack to deal with the fact
            % that x can have any number of dimensions.)
            extraDims = ndims(x) - 2;
            extraColons = repmat(':', 1, extraDims);
            extraColons = mat2cell(extraColons, 1, ones(1, extraDims));
            out = out(:, outCols(:), extraColons{:});

        elseif ( isa(x, 'chebfun') )
            % Call COMPOSE():
            out = compose(x, f);
            
        elseif ( isequal(x, ':') )
            % Return f:
            if ( (numel(columnIndex) == size(f, 2)) && ...
                 all(columnIndex == 1:size(f, 2)) )
                out = f;
            else
                % Extract the required columns:
                out = extractColumns(f, columnIndex);
            end
            
        else
            error('CHEBFUN:CHEBFUN:subsref:nonnumeric',...
              'Cannot evaluate chebfun for non-numeric type.')
          
        end
        
        % Deal with row CHEBFUN objects:
        if ( isTransposed )
            out = permute(out, [2, 1, 3:ndims(out)]);
        end
        
        % Recurse on SUBSREF():
        if ( numel(index) > 1 )
            out = subsref(out, index(2:end));
        end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(f, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '{}'

        if ( length(idx) == 1 )
            if ( isequal(idx{1}, ':') )
                % F{:} returns F:
                out = f;
            else
                error('CHEBFUN:CHEBFUN:subsref:badDomain', ...
                    'Invalid domain syntax.')
            end
            
        elseif ( size(idx, 1) == 1 )
            % F{s1,s2,...,sk} returns RESTRICT(F, [s1,s2,...,sk]):
            x = cat(2, idx{:});
            out = restrict(f, x);
            out = simplify(out);
            
        else
            error('CHEBFUN:CHEBFUN:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')
            
        end
        
    otherwise
        
        error('CHEBFUN:CHEBFUN:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end


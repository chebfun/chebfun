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
%   F{S1, S2} restricts F to the domain [S1, S2] < [F.ENDS(1), F.ENDS(end)]. See
%   CHEBFUN/RESTRICT for further details. Note that F{[S1, S2]} is not supported
%   due to the behaviour of the MATLAB subsref() command.
%
% See also FEVAL, COMPOSE, GET, RESTRICT, SUBSREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'
        
        % Deal with row CHEBFUN objects:
        isTransposed = f.isTransposed;
        if ( isTransposed )
            f = f.';
            if ( length(idx) > 1 )
                idx(1:2) = idx([2,1]);
            end
        end
        
        % Where to evaluate:
        x = idx{1}; 
        
        % Initialise:
        columnIndex = 1:size(f, 2); % Column index for array-valued CHEBFUNs.
        varin = {};                 % Additional arguments.
        
        % Deal with additional arguments:
        if ( (length(idx) == 2) && ...
             any(strcmpi(idx{2}, {'left', 'right', '-', '+'})) )
            % f(x, 'left') or f(x, 'right'):
            varin = {idx(2)};
            
        elseif ( (length(idx) == 2) && ...
                 ((max(idx{2}) <= columnIndex(end)) || strcmp(idx{2}, ':')) )
            % f(x, m), for array-valued CHEBFUN objects:
            columnIndex = idx{2};         
            
        elseif ( length(idx) == 2 && strcmp(idx{2}, ':') )
            % This is OK.
            
        elseif ( length(idx) > 1 )
            error('CHEBFUN:subsref:dimensions', ...
                'Index exceeds CHEBFUN dimensions.')
            
        end

        % Compute the output:
        if ( isnumeric(x) )
            % Call FEVAL():
            out = feval(f, x, varin{:});

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
            error('CHEBFUN:subsref:nonnumeric',...
              'Cannot evaluate chebfun for non-numeric type.')
          
        end
        
        % Deal with row CHEBFUN objects:
        % (NB:  We call PERMUTE instead of TRANSPOSE in case OUT is
        % multidimensional).
        if ( isTransposed )
            out = permute(out, [2 1 3:ndims(out)]);
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
                error('CHEBFUN:subsref:baddomain', 'Invalid domain syntax.')
            end
            
        elseif ( size(idx, 1) == 1 )
            % F{s1,s2,...,sk} returns RESTRICT(F, [s1,s2,...,sk]):
            x = cat(2, idx{:});
            out = restrict(f, x);
            
        else
            error('CHEBFUN:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')
            
        end
        
    otherwise
        
        error('CHEBFUN:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end


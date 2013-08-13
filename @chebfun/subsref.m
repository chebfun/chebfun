function varargout = subsref(f, index)
%SUBSREF   Chebfun subsref.
% ( )
%   F(X) returns the values of the chebfun F evaluated on the array X.
%
%   If X falls on a breakpoint of F, the corresponding value from F.IMPULSES is
%   returned. F(X, 'left') or F(X, 'right') will evaluate F to the left or right
%   of the breakpoint respectively.
%
%   F(G), where G is also a chebfun, computes the composition of F and G.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2} restricts F to the domain [S1, S2] < [F.ENDS(1), F.ENDS(end)].
%
% See also FEVAL, COMPOSE, GET, RESTRICT, SUBSREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Should we allow specific access to columns in array-valued CHEBFUN
% objects? For example, F(0, 2) to evaluate the 2nd columns at x = 0?

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        % Where to evaluate:
        x = idx{1}; 

        % Initialise:
        columnIndex = 1:size(f.funs{1}, 2); % For array-valued CHEBFUN objects.
        varin = {};                         % Additional arguments.
        
        % Deal with additional arguments:
        if ( length(idx) == 2 ) && ...
                ( any(strcmpi(idx{2}, {'left', 'right', '-', '+'})) )
            % f(x, 'left') or f(x, 'right'):
            varin = {idx(2)};
            
        elseif ( length(idx) == 2 ) && max(idx{2}) <= columnIndex(end)
            % f(x, m), for array-valued CHEBFUN objects:
            columnIndex = idx{2};         
            
        elseif ( length(idx) > 1 )
            error('CHEBFUN:subsref:dimensions',...
                'Index exceeds chebfun dimensions.')
            
        end

        % Compute the output:
        if ( isnumeric(x) )
            % Call FEVAL():
            out = feval(f, x, varin{:});
            out = out(:,columnIndex);
            
        elseif ( isa(x, 'chebfun') || isa(x, 'function_handle') )
            % Call COMPOSE():
            out = compose(f, x);
            
        elseif ( isequal(x, ':') )
            % Return f:
            out = f;
            
        else
            error('CHEBFUN:subsref:nonnumeric',...
              'Cannot evaluate chebfun for non-numeric type.')
          
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
                return
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
        
        error('CHEBFUN:UnexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end


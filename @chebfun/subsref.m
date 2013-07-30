function varargout = subsref(f, index)
% SUBSREF   Chebfun subsref.
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
%   F{S} restricts F to the domain [S(1) S(end)] < [F.ENDS(1) F.ENDS(end)].
%
% See also CHEBFUN/FEVAL, CHEBFUN/GET, CHEBFUN/RESTRICT

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    case '.'

        % Call get for .PROP access.
        if ( strcmp(idx, 'domain') && numel(index) > 1 )
            varargout{1} = f.domain(index(2).subs{:});
        elseif strcmp(idx, 'funs')
            varargout{1} = get(f, idx);
            if ( numel(index) > 1 )
                index(1) = [];
                varargout{1} = subsref(varargout{1}, index);
            end
        else
            varargout = get(f, idx);
            if ( ~iscell(varargout) )
                varargout = {varargout};
            end
        end

    case '()'

        x = idx{1}; % Where to evaluate.

        % Deal with extra arguments:
        if ( length(idx) == 1 )
            varin = {};
        elseif ( length(idx) == 2 )
            % f(x, 'left') or f(x, 'right')
            if ( any(strcmpi(idx{2}, {'left', 'right', '-', '+'})) )
                varin = {idx(2)};
            end
        elseif ( length(idx) ~= 1 )
            error('CHEBFUN:subsref:dimensions',...
                'Index exceeds chebfun dimensions.')
        end

        if ( isnumeric(x) )
            varargout = { feval(f, x, varin{:}) };
        elseif ( isa(x, 'chebfun') || isa(x, 'function_handle') )
            varargout = { compose(f, x) };
        elseif ( isequal(x, ':') )
            varargout = { f };
        else
            error('CHEBFUN:subsref:nonnumeric',...
              'Cannot evaluate chebfun for non-numeric type.')
        end

    case '{}'

        if ( length(idx) == 1 )
            if ( isequal(idx{1}, ':') )
                varargout = {f};
                return
            else
                error('CHEBFUN:subsref:baddomain',...
                    'Invalid domain syntax.')
            end
        elseif ( length(idx) == 2 )
            x = cat(2, idx{:});
        else
            error('CHEBFUN:subsasgn:dimensions',...
                'Index exceeds chebfun dimensions.')
        end
        varargout = { restrict(f, x) };
    otherwise
        error('CHEBFUN:UnexpectedType',...
            ['??? Unexpected index.type of ' index(1).type]);
end

% if ( length(index) > 1 )
%     varargout = {subsref([varargout{:}], index(2:end))};
% end

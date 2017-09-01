function f = compose(f, op, g, data, pref)
%COMPOSE   Composition of SINGFUN objects.
%   COMPOSE(F, OP) returns a ONEFUN representing OP(F) where F is a SINGFUN 
%   object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a ONEFUN representing OP(F, G) where at least one 
%   of F and G are SINGFUN objects, and OP is a function handle.
%
%   COMPOSE(F, G) returns a ONEFUN representing G(F), where at least of F and G 
%   are SINGFUN. If the range of F is not in [-1, 1] then an error is thrown.
%
%   COMPOSE(F, OP, G, DATA, PREF) or COMPOSE(F, OP, [], DATA, PREF) uses the
%   constructor data in the structure DATA and the options passed by the
%   CHEBFUNPREF or preference structure PREF to build the returned ONEFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if ( nargin < 5 )
    pref = chebfunpref();
else
    pref = chebfunpref(pref);
end

if ( nargin < 4 )
    data = struct();
end

if ( (nargin < 3) || isempty(g) )
    nfuns = 1;
else
    nfuns = 2;
end

if ( nfuns == 2 )
    if ( size(f, 2) ~= size(g, 2) )
        error('CHEBFUN:SINGFUN:compose:dim', ...
              'Matrix dimensions must agree.')
    end
    
elseif ( isa(op, 'smoothfun') || isa(op, 'singfun') )
    % If OP is a SMOOTHFUN or a SINGFUN, we check the dimension and range of F
    if ( (size(op, 2) > 1) && (size(f, 2) > 1) )
        error('CHEBFUN:SINGFUN:compose:arrval', ...
              'Cannot compose two array-valued ONEFUN objects.')
    end

    % [TODO]: How to check domain-range matchup in SINGFUN?  ('values' may not
    % be available for every tech.)
    % if ( norm(get(f, 'values'), inf) > 1 )
    %    error('CHEBFUN:SINGFUN:compose:range', ...
    %          'The range of f is not contained in the domain of g.')
    % end

end

if ( nfuns == 1 ) % OP(F)
    op = @(x) feval(op, feval(f, x));
else % OP(F, G)
    op = @(x) feval(op, feval(f, x), feval(g, x));
end

% Call SINGFUN constructor:
f = singfun(op, data, pref);

% Simplify:
f = simplify(f);

% Degenerate to SMOOTHFUN:
f = singFun2SmoothFun(f);

end

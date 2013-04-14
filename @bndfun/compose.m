function f = compose(f, op, g, pref)
%COMPOSE  Compostition of BNDFUN objects.
%   COMPOSE(F, OP) returns a BNDFUN representing OP(F) where F is also a
%   BNDFUN2 object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns OP(F, G) where F and G are BNDFUN objects,
%   and OP is a function handle.
%
%   COMPOSE(F, G) returns a BNDFUN representing G(F), where both F and G are
%   also BNDFUN objects. If the range of F is not in [-1, 1] then an error is
%   thrown.
%
%   COMPOSE(F, OP, PREF), COMPOSE(F, OP, G, PREF), or COMPOSE(F, G, PREF) uses
%   the options passed by the prefences structure PREF. In particular, one can
%   pass a PREF.(class(F)).refinmentFunction which takes advantage of the fact
%   that F (and possibly OP or G) are BNDFUN objects.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
if nargin == 2
    f.onefun = compose(f.onefun,op);
elseif nargin == 3
    if ( isstruct(g) )
        % Third argument passed was a preference structure.
        
        f.onefun = compose(f.onefun,op,g);
    else
        % Third argument passed was a bndfun
        if ( ~checkDomain(f, g))
            
            % Throw an error if domains don't match:
            error('BNDFUN:compose:domainMismatch',...
                'Bndfuns domains must agree.')
        end
        % Compose the onefun of f with the onefun of g:
        f.onefun = compose(f.onefun,op,g.onefun);
    end
else
    f.onefun = compose(f.onefun,op.g.onefun,pref);
end

end

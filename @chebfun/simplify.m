function F = simplify(F, tol)
%SIMPLIFY  Simplify a CHEBFUN.
%  G = SIMPLIFY(F) attempts to compute a CHEBFUN G which is a 'simplified'
%  version of F in that length(G) <= length(F), but ||G - F|| is small in a
%  relative sense: ||G - F|| < EPSLEVEL(G)*VSCALE(G).  The relative error
%  threshold tolerance is chosen based on F's own global accuracy estimate (via
%  F.VSCALE and F.EPSLEVEL) and the local VSCALEs of F's individual FUN
%  objects.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of the
%  default simplification tolerances as the relative threshold level.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

if ( nargin == 1 )
    tol = chebpref().techPrefs.eps;
elseif ( isa(tol, 'chebpref') )
    tol = tol.techPrefs.eps;
end

% Loop over the columns:
for k = 1:numel(F)
    F(k) = columnSimplify(F(k), tol);
end

end

function f = columnSimplify(f, tol)

% Choose a tolerance:
% [TODO]: This seems to be the best we can do without vector epslevels?
glob_acc = epslevel(f).*vscale(f); % Global error estimate of the CHEBFUN.
loc_vscl = get(f, 'vscale-local'); % Local vscale of the FUN objects.
tol = max(max(tol, glob_acc) ./ loc_vscl, [], 2); % Simplification tolerances.

% Simplfy each of the FUN objects:
for k = 1:numel(f.funs)
    f.funs{k} = simplify(f.funs{k}, tol(k,:));
end

end

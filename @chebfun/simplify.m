function f = simplify(f, tol)
%SIMPLIFY  Simplify a CHEBFUN.
%  G = SIMPLIFY(F) attempts to compute a 'simplified' version G such that
%  length(G) <= length(F), but ||G - F|| is small in a relative sense: ||G - F||
%  < epslevel(G)*vscale(G).
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses TOL instead of the
%  default CHEBFUN.PREF('EPS') preference as the relative threshold level for
%  deciding whether a coefficient is small enough to be zeroed. Here,
%  epslevel(G) is the maximum of epslevel(F) and TOL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin == 1 )
    tol = chebfun.pref('eps');
end

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

function F = simplify(F, tol)
%SIMPLIFY  Simplify a CHEBFUN.
%  G = SIMPLIFY(F) attempts to compute a CHEBFUN G which is a 'simplified'
%  version of F in that length(G) <= length(F), but ||G - F|| is small in a
%  relative sense: ||G - F|| < EPS*VSCALE(G).  The relative error threshold
%  tolerance is chosen based on F's own global accuracy estimate (via VSCALE(F)
%  and EPS) and the local VSCALEs of F's individual FUN objects.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses the scalar TOL instead
%  of the default simplification tolerance as the relative threshold level.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    tol = chebfunpref().techPrefs.eps;
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

% Loop over the columns:
for j = 1:numel(F)
    % Adjust tolerance based on ratio of the global and local vscales:
    vscaleLocal = get(F(j), 'vscale-local');
    vscaleGlobal = vscale(F(j));
    tolj = tol.*vscaleLocal./vscaleGlobal;

    for k = 1:numel(F(j).funs)
        F(j).funs{k} = simplify(F(j).funs{k}, tolj(k));
    end
end

end

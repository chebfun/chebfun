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

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    tol = [];
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

% Loop over the columns:
for k = 1:numel(F)
    F(k) = columnSimplify(F(k), tol);
end

end

function f = columnSimplify(f, tol)

if ( nargin < 2 || isempty(tol) )
    p = chebfunpref;
    tol = p.eps;
end

if ( length(tol) ~= numel(f.funs) )
    tol = repmat(tol, 1, numel(f.funs));
end

vscl = f.vscale;

% Simplfy each of the FUN objects:
for k = 1:numel(f.funs)
    scl = vscl./get(f.funs{k},'vscale');
    f.funs{k} = simplify(f.funs{k}, tol(k)*scl);
end

end

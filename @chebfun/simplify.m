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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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
    % Choose a tolerance:
    loc_epsl = get(f, 'epslevel-local'); % Epslevel of each columns and each fun
    loc_vscl = get(f, 'vscale-local');   % Vscale of each columns and each fun
    loc_acc = loc_epsl.*loc_vscl;        % Pointwise multiply of the two
    col_vscl = max(loc_vscl, [], 1);     % Vscale of each column
    col_vscl(col_vscl == 0) = 1;         % Remove zero vscales (TODO: Improve?)
    tol = bsxfun(@rdivide, loc_acc, col_vscl); % Factor out col_vscl
elseif ( isscalar(tol) )
    tol = repmat(tol, 1, numel(f.funs));
end

% Simplfy each of the FUN objects:
for k = 1:numel(f.funs)
    f.funs{k} = simplify(f.funs{k}, tol(k));
end

end

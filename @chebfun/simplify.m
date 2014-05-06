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
% See http://www.chebfun.org for Chebfun information.

if ( nargin == 1 )
    tol = chebfunpref().techPrefs.eps;
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

% Loop over the columns:
for k = 1:numel(F)
    F(k) = columnSimplify(F(k), tol);
end

end

function f = columnSimplify(f, tol)

if ( nargin < 2 ) 
    % Choose a tolerance:
    loc_epsl = get(f, 'epslevel-local'); % Epslevel of each columns and each fun
    loc_vscl = get(f, 'vscale-local');   % Vscale of each columns and each fun
    loc_acc = loc_epsl.*loc_vscl;        % Pointwise multiply of the two
    glob_acc = max(max(loc_acc));        % Maximum over all columns and funs
    col_vscl = max(loc_vscl, [], 1);     % Vscale of each column
    col_vscl(col_vscl == 0) = 1;         % Remove zero vscales (TODO: Improve?)
    tol = loc_acc./col_vscl;             % Factor out col_vscl (TODO: loc_vscl?)
end

% Simplfy each of the FUN objects:
for k = 1:numel(f.funs)
    f.funs{k} = simplify(f.funs{k}, tol);
end

end

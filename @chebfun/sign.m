function f = sign(f, varargin)
%SIGN    Sign function of a CHEBFUN.
%   G = SIGN(F) returns a piecewise constant CHEBFUN G such that G(x) = 1 in the
%   interval where F(x) > 0, G(x) = -1 in the interval where F(x) < 0 and G(x) =
%   0 in the interval where F(x) = 0. Breakpoints in G are introduced at zeros
%   of F.
%
%   For the nonzero elements of complex F, sign(F) = F.
%
% See also ABS, HEAVISIDE, ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Deal with the empty case:
if ( isempty(f) )
    return
end

% If f is not real, sign returns f?
if ( ~isreal(f) )
    % [TODO]: Perhaps it should return:
%     out = f./abs(f);
    return
end

% Set scales and tolerances:
vs = vscale(f);
hs = hscale(f);
el = epslevel(f);
tol = 100*el*vs;

oldDom = f.domain;
oldImps = f.impulses;
numFuns = numel(f.funs);
numCols = size(f.funs{1}, 2);
isTransposed = f.isTransposed;
f.isTransposed = false;

%% %%%%%%%%%%%%%%%%%%%%%%% DETERMINE BREAKPOINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the roots of the CHEBFUN:
r = roots(f, 'nozerofun', 'noimps');

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(r(:)) and remove any remaining NaNs.
r = unique(r(:));
r(isnan(r)) = [];

% Discard any roots that are closer than the accuracy of the CHEBFUN:
rtol = el*hs*vs;
r([ false ; (diff(r) < rtol) ]) = [];
% Also avoid introducing new breakpoints close to the ends of the domain:
r(any(abs(bsxfun(@minus, r, oldDom([1,end]))) < rtol, 2)) = [];

% Add the domain boundaries to the roots vector:
r = [f.domain(1) ; r ; f.domain(end)];

% Non-trivial impulses should _not_ be removed:
lvals = [get(f, 'lval-local') ; get(f, 'rval')];
rvals = [get(f, 'lval') ; get(f, 'rval-local')];
% Check that the sign of the values is different _and_ values are far apart:
idxl = sign(oldImps(:,:,1)) ~= sign(lvals) & abs(oldImps(:,:,1) - lvals) > 100*tol;
idxr = sign(oldImps(:,:,1)) ~= sign(rvals) & abs(oldImps(:,:,1) - rvals) > 100*tol;
% Look for a jump to left or right in any column: 
idx = any(idxl | idxr, 2); % idx is st exists a non-trivial jump at oldDom(idx).
% Append the new breakpoints with nontrival impulses to r and sort:
[r2, ignored, idx2] = unique([oldDom(idx).' ; r]);
% Location of existing breakpoints inside new r vector:
idx2 = idx2(1:sum(idx));

%% %%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCT NEW CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate at an arbitrary point between each root:
c = 0.5912;
nr = length(r);
x = c*r(1:nr-1) + (1-c)*r(2:nr);
vals = sign( feval(f, x(:)) );

% Build new CHEBFUN:
f.domain = r.';
f.funs = cell(1, nr-1);
for k = 1:nr-1
    f.funs{k} = fun.constructor(vals(k,:), f.domain(k:k+1));
end

% Reassign impulses:
f.impulses = chebfun.jumpVals(f.funs);
f.impulses(idx2,:) = sign(oldImps(idx,:,1));

% Merge g to simplify unecessary breakpoints:
f = merge(f, find(idx2));

f.isTransposed = isTransposed;

end

function H = hypot(f, g, pref)
%HYPOT   Robust computation of the square root of the sum of squares.
%   H = HYPOT(F, G) returns SQRT(ABS(F).^2 + ABS(G).^2) for two CHEBFUN objects
%   F and G (or a CHEBFUN and a double) carefully computed to avoid underflow
%   and overflow.
%
% Example:
%       f = chebfun(@(x) 3*[1e300*x 1e-300*x]);
%       g = chebfun(@(x) 4*[1e300*x 1e-300*x]);
%       % h1 = sqrt(f.^2 + g.^2) % This will fail because of overflow
%       h2 = hypot(f, g)
%
% See also ABS, NORM, SQRT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Grab some preferences:
if ( nargin < 3 )
    pref = chebfun.pref();
end

% Insert breaks at the roots:
if ( isa(f, 'chebfun') )
    f = breaksAtRoots(f);
end
if ( isa(g, 'chebfun') )
    g = breaksAtRoots(g);
end

% Call compose:
H = compose(f, @hypot, g, pref);

end

function f = breaksAtRoots(f)
% Insert breaks whre the function f is zero (with some tolerance).

% [TODO]: This is used in other places (such as ABS) and needs de-duplicating.

% Locate roots:
r = roots(f, 'nozerofun');

% Since each column of an array-valued CHEBFUN must have the same breakpoints,
% we simply take unique(r(:)) and remove any remaining NaNs.
r = unique(r(:));
r(isnan(r)) = [];

% Discard any roots which are closer than the accuracy of the CHEBFUN:
el = epslevel(f);
vs = min(vscale(f),1);
hs = hscale(f);
rtol1 = el*hs*vs;
r([false ; diff(r) < rtol1]) = [];

% Avoid introducing new breakpoints close to an existing ones:
rtol2 = el*vs*max(min(diff(f.domain)), 1);
r(any(abs(bsxfun(@minus, r, f.domain)) < rtol2, 2)) = [];

% Add new breaks if required:
if ( ~isempty(r) )
    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = unique([f.domain, r.']);

    % Introduce these breakpoints to f:
    f = restrict(f, dom);
end

end
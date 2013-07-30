function H = hypot(f, g, pref)
%HYPOT  Robust computation of the square root of the sum of squares.
%   H = HYPOT(F, G) returns SQRT(ABS(F).^2 + ABS(G).^2) for two chebfun objects
%   F and G (or a chebfun and a double) carefully computed to avoid underflow
%   and overflow.
%
% Example:
%        x = chebfun('x', [-1, 1]);
%        f = 3*[1e300*x 1e-300*x];
%        g = 4*[1e300*x 1e-300*x];
%        % h1 = sqrt(f.^2 + g.^2) % This will fail because of overflow
%        h2 = hypot(f, g)
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
    f = rootsAtBreaks(f, pref);
end
if ( isa(g, 'chebfun') )
    g = rootsAtBreaks(g, pref);
end

% Call compose:
H = compose(f, @hypot, g, pref);

end

function f = rootsAtBreaks(f, pref)
% Insert breaks whre the function f is zero (with some tolerance).

% Abs is singular at roots, so locate these:
r = roots(f, 'nozerofun');
% Avoid adding new breaks where not needed:
if ( ~isempty(r) )

    % Choose a tolerance:
    rtol = 100*pref.chebfun.eps.*max(min(diff(f.domain)), 1);

    % Remove if sufficiently close to an existing break points:
    [rem, ignored] = find(abs(bsxfun(@minus, r , f.domain)) < rtol); %#ok<NASGU>
    r(rem) = [];

end

% Add new breaks if required:
if ( ~isempty(r) )
    % Get the domain with the new breakpoints: (union is not required, by above)
    dom = sort([f.domain, r']);

    % Introduce these breakpoints to f:
    f = restrict(f, dom);
end

end
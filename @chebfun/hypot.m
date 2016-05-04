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

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab some preferences:
if ( nargin < 3 )
    pref = chebfunpref();
end

% Insert breaks at the roots:
if ( isa(f, 'chebfun') )
    f = addBreaksAtRoots(f, pref);
end
if ( isa(g, 'chebfun') )
    g = addBreaksAtRoots(g, pref);
end

% Call compose:
H = compose(f, @hypot, g, pref);

end

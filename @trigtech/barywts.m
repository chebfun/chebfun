function v = barywts(n)
%BARYWTS  Barycentric weights for trigpts.
%   BARYWTS(N) returns the N barycentric weights for trigonometric-polynomial
%   interpolation on an equispaced grid [1].
%
% See also TRIGTECH.BARY, TRIGPTS. 

% Reference:
%
% [1] P. Henrici, Barycentric formulas for interpolating trigonometric
% polynomials and their conjugates, Numer. Math., 33 (1979), pp. 225–234.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )                      % Special case (no points)
    v = [];
elseif ( n == 1 )                  % Special case (single point)
    v = 1;
else                               % General case.
    v = ones(n, 1);
    v(2:2:end) = -1;
end

end

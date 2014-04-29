function v = barywts(n)
%BARYWTS  Barycentric weights for Chebyshev points of 1st kind.
%   BARYWTS(N) returns the N barycentric weights for polynomial interpolation on
%   a Chebyshev grid of the 1st kind.
%
% See also BARY, CHEBPTS. 

% Reference:
%
% [1] Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange
% Interpolation", SIAM Review 46 (2004), pp. 501-517.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )                      % Special case (no points)
    v = [];
elseif ( n == 1 )                  % Special case (single point)
    v = 1;
else
    v = sin((2*(n-1:-1:0)+1)*pi/(2*n)).';
    v(end-1:-2:1) = -v(end-1:-2:1);
end

end

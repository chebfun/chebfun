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

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( n == 0 )                      % Special case (no points)
    v = [];
elseif ( n == 1 )                  % Special case (single point)
    v = 1;
else
    v = sin(((n-1:-1:0)+.5)*pi/n).';
    % The following flipping trick forces symmetry. Also due to the nature of 
    % the sine function, those computed with a big argument are replaced by ones
    % with a small argument, improving the relative accuracy.
    v(1:floor(n/2)) = v(end:-1:n-floor(n/2)+1);
    v(end-1:-2:1) = -v(end-1:-2:1);
end

end

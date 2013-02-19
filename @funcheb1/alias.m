function coeffs = alias(coeffs, m)
%ALIAS  Alias Chebyshev coefficients of the 1st kind.
%   C = ALIAS(C, M) aliases the Chebyshev coefficients C to have length M. If
%   length(C) > M, this is equivalent to padding with zeros.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Reference:
% 
% [1] L. N. Trefethen, Approximation Theory and Approximation Practice,
% SIAM, 2013.

n = size(coeffs, 1);

% Pad with zeros:
if ( m > n )
    coeffs = [zeros(m-n, size(coeffs, 2)) ; coeffs];
    return
end

% Actually alias (eq. (4.4) of [1]):
coeffs = coeffs(end:-1:1,:);
for j = (m + 1):n
    k = abs( mod( j + m - 3 , 2*m - 2 ) - m + 2 ) + 1;
    coeffs(k, :) = coeffs(k, :) + coeffs(j, :);
end
coeffs = coeffs(m:-1:1, :);

end
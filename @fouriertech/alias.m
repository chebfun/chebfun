function coeffs = alias(coeffs, m)
%ALIAS   Alias Fourier coefficients on equally spaced grid.
%   ALIAS(C, M) aliases the Fourier coefficients stored in the column vector C
%   to have length M. If M > LENGTH(C), the coefficients are padded with zeros.
%   If C is a matrix of coefficients, each of the columns is aliased to length
%   M.
%
% See also PROLONG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful References:  TODO: update these
%
%   L.N. Trefethen, Approximation Theory and Approximation Practice, SIAM, 2013
%   Page 27.
%
%   Fox, L. and Parker, I. B., Chebyshev polynomials in Numerical Analysis,
%   Oxford University Press, 1972.  (pp. 67)
%
%   Mason, J. C. and Handscomb, D. C., Chebyshev polynomials, Chapman &
%   Hall/CRC, Boca Raton, FL, 2003.  (pp. 153)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(coeffs, 1);

% Pad with zeros:
if ( m > n )
%     coeffs = [ zeros(m-n, size(coeffs, 2)) ; coeffs ];
    z = zeros(ceil((m-n)/2), size(coeffs, 2));
    coeffs = [ z ; coeffs ; z ];
    return
end

% TODO: We're just going to truncate, rather than alias.
n = size(coeffs, 1);
n2 = floor(n/2);
coeffs(end-n2+1+m:end, :) = [];
coeffs(1:n2-m, :) = [];

end

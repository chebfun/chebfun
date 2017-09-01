function coeffs = alias(coeffs, m)
%ALIAS   Alias Chebyshev coefficients on the 2nd kind Chebyshev grid.
%   ALIAS(C, M) aliases the Chebyshev coefficients stored in the column vector C
%   to have length M. If M > LENGTH(C), the coefficients are padded with zeros.
%   If C is a matrix of coefficients, each of the columns is aliased to length
%   M.
%
% See also PROLONG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful References:
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
    coeffs = [ coeffs; zeros(m-n, size(coeffs, 2)) ];
    return
end

% Alias coefficients: (see eq. (4.4) of Trefethen, Approximation Theory and
% Approximation Practice, SIAM, 2013):
if ( m == 1 )
    % Reduce to a single point:
    e = ones(1, ceil(n/2)); 
    e(2:2:end) = -1;
    coeffs = e*coeffs(1:2:end,:);
elseif ( m > n/2 )
    % If m > n/2, only single coefficients are aliased, and we can vectorise.
    j = (m + 1):n;
    k = abs(mod(j + m - 3, 2*m - 2) - m + 2) + 1;
    coeffs(k,:) = coeffs(k,:) + coeffs(j,:);
else
    % Otherwise we must do everything in a tight loop. (Which is slower!)
    for j = (m + 1):n
        k = abs(mod(j + m - 3, 2*m - 2) - m + 2) + 1;        
        coeffs(k,:) = coeffs(k,:) + coeffs(j,:);
    end
end

% Truncate:
coeffs = coeffs(1:m,:);

end

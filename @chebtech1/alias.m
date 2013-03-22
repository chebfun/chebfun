function coeffs = alias(coeffs, m)
%ALIAS   Alias Chebyshev coefficients on the 1st kind Chebyshev grid.
%   ALIAS(C, M) aliases the Chebyshev coefficients stored in the column
%   vector C to have length M. If M > LENGTH(C), the coefficients are 
%   padded with zeros. If C is a matrix of coefficients, each of the 
%   columns is aliased to length M.

% Note that aliasing for the 1st kind Chebyshev grid is different from its
% counterpart for 2nd kind Chebyshev grid, though they both live in the 
% coefficient space with respect to Chebyshev polynomials of 1st kind. 

%   References:
%
%   [1] Fox, L. and Parker, I. B., Chebyshev polynomials in Numerical 
%   Analysis, Oxford University Press, 1972.  (pp. 67)
%
%   [2] Mason, J. C. and Handscomb, D. C., Chebyshev polynomials, Chapman 
%   & Hall/CRC, Boca Raton, FL, 2003.  (pp. 153)
%
%   [3] Chebfun 5v Working Note, 03/05/2013.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = size(coeffs, 1);
n2 = size(coeffs, 2);

% Pad with zeros:
if ( m > n )
    coeffs = [ zeros(m-n, size(coeffs, 2)); coeffs ];
    return
end

% It's more natural to work with the coefficients in the other order:
coeffs = coeffs(end:-1:1,:);

% Alias coefficients (see discussion above):
if ( m == 1 )
    % Reduce to a single point:
    e = ones(1, ceil(n/2)); 
    e(2:2:end) = -1;
    coeffs = e*coeffs(1:2:end,:);
elseif ( m > n/2 )
    % If m > n/2, only single coefficients are aliased, and we can vectorise.
    j = ((m + 1):n).';
    k = abs(mod(j + m - 2, 2*m) - m + 1) + 1;
    p = floor((j-1+m)/(2*m));
    t = (-1).^p;
    coeffs(k,:) = coeffs(k,:) + repmat(t,1,n2).*coeffs(j,:);
else
    % Otherwise we must do everything in a tight loop. (Which is slower!)
    for j = (m + 1):n
        k = abs(mod(j + m - 2, 2*m) - m + 1) + 1;
        p = floor((j-1+m)/(2*m));
        t = (-1)^p;
        coeffs(k,:) = coeffs(k,:) + t*coeffs(j,:);
    end
end

% Flip the coefficients back again:
coeffs = coeffs(m:-1:1,:);

end

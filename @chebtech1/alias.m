function coeffs = alias(coeffs, m)
%ALIAS   Alias Chebyshev coefficients on the 1st kind Chebyshev grid.
%   ALIAS(C, M) aliases the Chebyshev coefficients stored in the column vector C
%   to have length M. If M > LENGTH(C), the coefficients are padded with zeros.
%   If C is a matrix of coefficients, each of the columns is aliased to length
%   M.
%
% See also PROLONG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the formula for aliasing on the 1st-kind Chebyshev grid is
% different from that for the 2nd-kind grid, even though the coefficients 
% being aliased are for 1st-kind Chebyshev polynomials in both cases. 
%
% Useful References:
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
    coeffs = [ coeffs ; zeros(m-n, size(coeffs, 2)) ];
    return
end

% Alias coefficients:
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
    coeffs(k,:) = coeffs(k,:) + bsxfun(@times, t, coeffs(j,:));
else
    % Otherwise we must do everything in a tight loop. (Which is slower!)
    for j = (m + 1):n
        k = abs(mod(j + m - 2, 2*m) - m + 1) + 1;
        sgn = 1 - 2*mod(floor((j - 1 + m)/(2*m)), 2);
        coeffs(k,:) = coeffs(k,:) + sgn*coeffs(j,:);
    end
end

% Truncate:
coeffs = coeffs(1:m,:);

end

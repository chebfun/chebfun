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
    k = ceil((m-n)/2);
    z = zeros(k, size(coeffs, 2));
    % Need to handle the odd vs. even case separately.
    if mod(n,2) == 0
        %
        % First account for the asymmetry in the coefficients when n is even.
        %
        % This will account for the cos(N/2) coefficient, which is stored
        % in the coeffs(n,:) entry, using properties of the complex
        % exponential.
        coeffs = [coeffs(n,:)/2;coeffs(1:n-1,:);coeffs(n,:)/2];
        coeffs = [z(1:end-1,:); coeffs; z];
        
        % Next check if m is odd.  If it is then coeffs is too long and we
        % need to remove the last row.
        if mod(m,2) == 1
            coeffs = coeffs(1:end-1,:);
        end
    else
        % There is no asymmetry in the coefficients just pad them.
        coeffs = [ z ; coeffs ; z];
        % Only need to check if m is even, in which case coeffs is too 
        % long and we need to remove the first row.
        if mod(m,2) == 0
            coeffs = coeffs(2:end,:);
        end
    end
%     numCoeffs = 2*k + n;
%     
%     % If the number of coefficients does not match m then remove the last
%     % row of the matrix.  This could happen when n is even and m is odd or
%     % when n is off and m is even.
%     if numCoeffs == m+1
%         coeffs = coeffs(2:end,:);
%     end
        
    return
end

% TODO: We're just going to truncate, rather than alias.
n = size(coeffs, 1);
n2 = floor(n/2);
coeffs(end-n2+1+m:end, :) = [];
coeffs(1:n2-m, :) = [];

end

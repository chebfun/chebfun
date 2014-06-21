function coeffs = alias(coeffs, m)
%ALIAS   Alias Fourier coefficients on equally spaced grid.
%   ALIAS(C, M) aliases the Fourier coefficients stored in the column vector C
%   to have length M. If M > LENGTH(C), the coefficients are padded with zeros.
%   If C is a matrix of coefficients, each of the columns is aliased to length
%   M.
%
% See also PROLONG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = size(coeffs, 1);

% Pad with zeros:
if ( m > n )
    
    k = ceil((m-n)/2);
    z = zeros(k, size(coeffs, 2));
    
    % Need to handle the odd vs. even case separately.
    if mod(n, 2) == 0
        
        % First, account for the asymmetry in the coefficients when n is even.
        % This will account for the cos(N/2) coefficient, which is stored
        % in the coeffs(n,:) entry, using properties of the complex
        % exponential.
        coeffs = [coeffs(n,:)/2;coeffs(1:n-1,:);coeffs(n,:)/2];
        coeffs = [z(1:end-1,:); coeffs; z];
        
        % Next, check if m is odd.  If it is then coeffs is too long and we
        % need to remove the last row.
        if mod(m, 2) == 1
            coeffs = coeffs(1:end-1,:);
        end
        
    else
        
        % There is no asymmetry in the coefficients, just pad them.
        coeffs = [ z; coeffs; z ];
        % Only need to check if m is even, in which case coeffs is too 
        % long and we need to remove the first row.
        if ( mod(m, 2) == 0 )
            coeffs = coeffs(2:end,:);
        end
    end
    
    return
    
end

isOdd = mod(n,2);

% Extract coefficients and flip then to start from lower modes (in
% absolute value) to higher modes.
if ( isOdd )
    constCoeffs = coeffs((n-1)/2+1,:);
    posCoeffs = coeffs((n-1)/2:-1:1,:);
    negCoeffs = coeffs((n-1)/2+2:n,:);
    
    if ( m == 1 )
        % Reduce to a single point:
        e = ones(1, ceil((n-1)/2));
        e(1:2:end) = -1;
        coeffs = constCoeffs + (e*posCoeffs + e*negCoeffs);
    else
        m2 = (m-1)/2;
        n2 = (n-1)/2;
        % It's more natural to work with the coefficients in the other order:
        coeffs = coeffs(end:-1:1,:);
        newCoeffs = 0*coeffs(n2-m2+1:n2+m2+1,:);
%         newCoeffs(m2+1,:) = 2*newCoeffs(m2+1,:);
        for j = -n2:0
            k = mod(j+m2+1,-m) + m2;
            coeffIndexK = k+m2+1;
            ceoffIndexJ = j+n2+1;
            newCoeffs(coeffIndexK,:) = newCoeffs(coeffIndexK,:) - (-1)^(k*j)*coeffs(ceoffIndexJ,:);
            coeffIndexK = -k+m2+1;
            ceoffIndexJ = -j+n2+1;
            newCoeffs(coeffIndexK,:) = newCoeffs(coeffIndexK,:) - (-1)^(k*j)*coeffs(ceoffIndexJ,:);
        end
%         newCoeffs(m2+1,:) = 0.5*newCoeffs(m2+1,:);
        coeffs = flipud(newCoeffs);
    end
end
% 
% % Simple solution
% x = fourtech.fourpts(m);
% values = feval(f,x);
% 
% if ( m == 1 )
%     % Reduce to a single point:
%     e = ones(1, ceil(n/2)); 
%     e(2:2:end) = -1;
%     coeffs = e*coeffs(1:2:end,:);
% elseif ( m > n/2 )
%     % If m > n/2, only single coefficients are aliased, and we can vectorise.
%     j = (m + 1):n;
%     k = abs(mod(j + m - 3, 2*m - 2) - m + 2) + 1;
%     coeffs(k,:) = coeffs(k,:) + coeffs(j,:);
% else
%     % Otherwise we must do everything in a tight loop. (Which is slower!)
%     for j = (m + 1):n
%         k = abs(mod(j + m - 3, 2*m - 2) - m + 2) + 1;
%         coeffs(k,:) = coeffs(k,:) + coeffs(j,:);
%     end
% end
% 
% % For now we just chop off the unwanted coefficients.
% % [TODO]: Use the aliasing formula for Fourier coefficients.
% 
% % This code simply chops off the unwanted coefficients.
% n = size(coeffs, 1);
% n2 = floor(n/2);
% % Negative modes to remove:
% negModeIndex = (n + 1 - n2 + floor(m/2)):n;
% % Positive modes to remove
% posModeIndex = 1:(n2-ceil(m/2)+mod(n, 2));
% coeffs(negModeIndex,:) = [];
% coeffs(posModeIndex,:) = [];

end
function coeffs = alias(coeffs, m)
%ALIAS   Alias Fourier coefficients on equally spaced grid.
%   ALIAS(C, M) aliases the Fourier coefficients stored in the column vector C
%   to have length M. If M > LENGTH(C), the coefficients are padded with zeros.
%   If C is a matrix of coefficients, each of the columns is aliased to length
%   M.
%
% See also PROLONG.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the number of coefficients.
n = size(coeffs, 1);

% Pad with zeros:
if ( m > n )
    
    k = ceil((m-n)/2);
    z = zeros(k, size(coeffs, 2));
    
    % Need to handle the odd vs. even case separately.
    if ( ~mod(n, 2) ) % n even.
        
        % First, account for the asymmetry in the coefficients when n is even.
        % This will account for the cos(N/2) coefficient, which is stored
        % in the coeffs(n,:) entry, using properties of the complex
        % exponential.
        coeffs = [ coeffs(1,:)/2; coeffs(2:n,:); coeffs(1,:)/2 ];
        coeffs = [ z; coeffs; z(1:end-1,:) ];
        
        % Next, check if m is odd. If it is, then coeffs is too long and we
        % need to remove the first row (the lowest degree
        % coefficients).
        if ( mod(m, 2) ) % m odd.
            coeffs = coeffs(2:end,:);
        end
        
    else % n odd.
        
        % There is no asymmetry in the coefficients, just pad them.
        coeffs = [ z; coeffs; z ];
        
        % Only need to check if m is even, in which case coeffs is too 
        % long and we need to remove the last row (the highest degree
        % coefficients).
        if ( ~mod(m, 2) ) % m even.
            coeffs = coeffs(1:end-1,:);
        end
    end
    
    return
    
end

% If the number of coefficients is even then we extend them by one so they
% are odd by exploiting the symmetry property.  This makes the code below a
% little cleaner since fewer cases need to be handled.
if ( ~mod(n, 2) ) % n even.
    coeffs(1,:) = 0.5*coeffs(1,:);
    coeffs = [ coeffs; coeffs(1,:) ];
    n = n + 1;
end

% Need to handle the odd and even case (for m) differently.
if ( mod(m, 2) ) % m is odd.

    if ( m == 1 )

        % Reduce to a single point:
        constCoeffs = coeffs((n-1)/2+1,:);
        posCoeffs = coeffs((n-1)/2:-1:1,:);
        negCoeffs = coeffs((n-1)/2+2:n,:);
        e = ones(1, ceil((n-1)/2));
        e(1:2:end) = -1;
        coeffs = constCoeffs + (e*posCoeffs + e*negCoeffs);

    else
        
        m2 = (m-1)/2;
        n2 = (n-1)/2;
        % Extract coefficients:
        aliasedCoeffs = coeffs(n2-m2+1:n2+m2+1,:);
        
        % The code below aliases the coefficients from the higher modes
        % onto the lower modes. The principle behind the formula is figure
        % out which of the higher Fourier modes are indistinguishable from
        % the lower Fouirer modes on the grid consisting of m equally
        % spaced points from [-1,1). In general, when m and n are odd the
        % following will be equal for j=-(n-1)/2 to -(m+1)/2
        % exp(1i*pi*j*x) = sgn*exp(1i*pi*k*x) where
        % k = mod(j+(m+1)/2, -m) + (m-1)/2 and sgn = (-1)^mod(j+k, 2);
        for j = -n2:-m2-1
            k = mod(j+m2+1, -m) + m2;
            coeffIndexK = k + m2 + 1;  % Index into aliasedCoeffs for mode k.
            coeffIndexJ = j + n2 + 1;  % Index into coeffs for mode j.
            sgn = (-1)^mod(j+k, 2);
            aliasedCoeffs(coeffIndexK,:) = aliasedCoeffs(coeffIndexK,:) + sgn*coeffs(coeffIndexJ,:);
            coeffIndexK = -k + m2 + 1;
            coeffIndexJ = -j + n2 + 1;
            aliasedCoeffs(coeffIndexK,:) = aliasedCoeffs(coeffIndexK,:) + sgn*coeffs(coeffIndexJ,:);
        end
        coeffs = aliasedCoeffs;

    end
    
else % m is even.
    
    m2 = m/2;
    n2 = (n-1)/2;
    
    % Put the coefficient for the cos(m/2*pi*x) in the first entry of the
    % coefficient vector.  This corresopnds to the exp(-1i*pi*m/2*x).  
    aliasedCoeffs = coeffs(n2+1-m2:n2+m2,:);
    % Extend the aliased cofficient vector so it is symmetric.  This allows
    % us to easily account for aliasing on the exp(-1i*pi*m/2*x) and 
    % exp(1i*pi*m/2*x) term in the code below without the need for a 
    % special if check.  The negative allows contribution for sin(m/2*pi*x)
    % to be removed.    
    aliasedCoeffs = [ aliasedCoeffs; -aliasedCoeffs(1,:) ];
    
    % Follow a similar structure to the odd m case above.  The main
    % difference is that the higher order modes do not change sign when
    % aliased onto an even point grid.
    for j = -n2:-m2
        k = mod(j+m2,-m) + m2;
        coeffIndexK = k + m2 + 1;  % Index into aliasedCoeffs for mode k.
        coeffIndexJ = j + n2 + 1;  % Index into coeffs for mode j.
        aliasedCoeffs(coeffIndexK,:) = aliasedCoeffs(coeffIndexK,:) + coeffs(coeffIndexJ,:);
        coeffIndexK = -k + m2 + 1;
        coeffIndexJ = -j + n2 + 1;
        aliasedCoeffs(coeffIndexK,:) = aliasedCoeffs(coeffIndexK,:) + coeffs(coeffIndexJ,:);
    end
    % Collapse the aliased coefficient vector down to an even number of
    % terms by adding in the aliasing of the exp(1i*pi*m/2*x) terms to the
    % exp(-1i*pi*m/2*x).
    aliasedCoeffs(1,:) = aliasedCoeffs(1,:) + aliasedCoeffs(end,:);
    aliasedCoeffs(end,:) = [];

    coeffs = aliasedCoeffs;

end

end

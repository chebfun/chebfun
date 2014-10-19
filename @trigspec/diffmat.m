function D = diffmat(N, m)
%DIFFMAT   Differentiation matrices for TRIGSPEC.
%   D = DIFFMAT(N, M) returns the differentiation matrix that takes N
%   Fourier coefficients and returns N coefficients that represent the mth 
%   derivative of the Fourier series. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs.
if ( nargin == 1 )
    m = 1;
end

% Create the differentation matrix.
if ( m > 0 )
    if ( mod(N, 2) == 0) % N even
        if ( mod(m, 2) == 1 ) % m odd
            D = (1i)^m*diag([0, -N/2+1:1:N/2-1]).^m;
        else % m even
            D = (1i)^m*diag(-N/2:1:N/2-1).^m;
        end
    else % N odd
        D = (1i)^m*diag(-(N-1)/2:1:(N-1)/2).^m;
    end
else
    D = eye(N);
end

end

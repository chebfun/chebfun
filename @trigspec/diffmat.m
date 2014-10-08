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
    if ( mod(N, 2) == 0)
        D = (1i)^m*diag([-N/2:1:N/2-1]).^m;
    else
        D = (1i)^m*diag([-(N-1)/2:1:(N-1)/2]).^m;
    end
else
    D = eye(N);
end

end

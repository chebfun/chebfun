function D = diffmat(N, m, flag)
%DIFFMAT   Differentiation matrices for TRIGSPEC.
%   D = DIFFMAT(N, M) returns the differentiation matrix that takes N
%   Fourier coefficients and returns N coefficients that represent the mth 
%   derivative of the Fourier series. 
%
%   D = DIFFMAT(N, M, 1) returns the same as DIFFMAT(N, M) unless N is even
%   and M is odd. If N is even and M is odd it returns the same matrix as 
%   DIFFMAT(N, M) except the (1,1) entry is (-1i*N/2)^m instead of 0. This 
%   flag is currently only used in the spherefun.poisson and 
%   spherefun.helmholtz commands.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs.
if ( nargin == 1 )
    m = 1;
elseif ( nargin < 3 )
    flag = 0;
end

% Create the differentation matrix.
if ( m > 0 )
    if ( mod(N, 2) == 0 ) % N even
        if ( mod(m, 2) == 1 ) % m odd
            D = (1i)^m*spdiags([0, -N/2+1:1:N/2-1]', 0, N, N).^m;
            if ( flag ) 
                % Set the (1,1) entry to (-1i*N/2)^m, instead of 0:
                D(1,1) = (-1i*N/2)^m;
            end
        else % m even
            D = (1i)^m*spdiags((-N/2:1:N/2-1)', 0, N, N).^m;
        end
    else % N odd
        D = (1i)^m*spdiags((-(N-1)/2:1:(N-1)/2)', 0, N, N).^m;
    end
elseif ( m == 0 )
    D = speye(N);
end

end

function c = jaccoeffs(f, n, alp, bet)
%JACCOEFFS   Compute Legendre series coefficients of a CHEBTECH object.
%   B = JACCOEFFS(F) returns the Legendre series coefficients of CHEBTECH F, so
%   that F = B(1)*P_0 + ... + B(N)*P_(N-1), where P_k is the kth Legendre
%   polynomial. B is a vector of the same length as that of F.
%
%   B = JACCOEFFS(F, N) returns the first N coefficients. If length(F) < N
%   then the additional entries of B are padded with zeros.
%
%   If F is an array-valued CHEBTECH, then a matrix of coefficients is returned
%   so that F(:,k) = B(1,k)*P_0 + ... + B(N,k)*P_(N-1).
%
% See also CHEBCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 3 )
    % Choose n = length(f)
    bet = alp;
    alp = n;
    n = length(f);
end

% Convert from Chebyshev to Jacobi coefficients.
alp
c = cheb2jac(f.coeffs, alp, bet);

% Truncate / pad to length n:
if ( nargin > 1 )
    s = size(c);
    if ( s(1) > n )
        c = c(1:n, :);
    else
        c = [c ; zeros(n-s(1), s(2))];
    end
end

end
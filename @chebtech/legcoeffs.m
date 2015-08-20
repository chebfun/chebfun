function b = legcoeffs(f, n)
%LEGCOEFFS   Compute Legendre series coefficients of a CHEBTECH object.
%   B = LEGCOEFFS(F) returns the Legendre series coefficients of CHEBTECH F, so
%   that F = B(1)*P_0 + ... + B(N)*P_(N-1), where P_k is the kth Legendre
%   polynomial. B is a vector of the same length as that of F.
%
%   B = LEGCOEFFS(F, N) returns the first N coefficients. If length(F) < N
%   then the additional entries of B are padded with zeros.
%
%   If F is an array-valued CHEBTECH, then a matrix of coefficients is returned
%   so that F(:,k) = B(1,k)*P_0 + ... + B(N,k)*P_(N-1).
%
% See also CHEBCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = cheb2leg(f.coeffs);

if ( nargin > 1 )
    s = size(b);
    if ( s(1) > n )
        b = b(1:n, :);
    else
        b = [b; zeros(n-s(1), s(2))];
    end
end

end

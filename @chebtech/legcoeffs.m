function b = legcoeffs(f, n)
%LEGCOEFFS   Compute Legendre series coefficients of a CHEBTECH object.
%   B = LEGCOEFFS(F) returns the Legendre series coefficients of CHEBTECH F, so
%   that F = B(N+1)*P_N + ... + B(1)*P_0, where P_k is the kth Legendre
%   polynomial.
%
%   B = LEGCOEFFS(F, N) returns the first N+1 coefficients. If length(F) < N + 1
%   then the additional entries of B are padded with zeros.
%
%   If F is an array-valued CHEBTECH, then a matrix of coefficients is returned
%   so that F(:,k) = B(N+1,k)*P_N + ... + B(1,k)*P_0.
%
% See also CHEBCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

b = cheb2leg(f.coeffs);

if ( nargin > 1 )
    s = size(b);
    if ( s(1) > n + 1 )
        b = b(end-n+1:end, :);
    else
        b = [zeros(n-s(1), s(2)) ; b];
    end
end

end


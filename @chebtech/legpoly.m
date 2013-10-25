function b = legpoly(f)
%LEGPOLY   Compute Legendre series coefficients of a CHEBTECH object.
%   B = LEGPOLY(F) returns the Legendre series coefficients of CHEBTECH F, so
%   that F = B(N+1)*P_N + ... + B(1)*P_0, where P_k is the kth Legendre
%   polynomial.
%
%   If F is an array-valued CHEBTECH, then a matrix of coefficients is returned
%   so that F(:,k) = B(N+1,k)*P_N + ... + B(1,k)*P_0.
%
% See also CHEBPOLY.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

c_cheb = f.coeffs;
b = zeros(size(c_cheb));
for k = 1:size(c_cheb, 2)
    b(:,k) = chebtech2.cheb2leg(c_cheb(:,k));
end
b = b.';

end


function c = jaccoeffs(f, varargin)
%LEGCOEFFS   Compute Legendre series coefficients of a BNDFUN object.
%   C = JACCOEFFS(F, N, A, B) returns the Legendre series coefficients of BNDFUN
%   F, so that F = C(N+1)*P^{(A,B)}_N + ... + C(1)*P^{(A,B)}_0, where
%   P^{(A,B)}_k is the degree k Jacobi polynomial (scaled to the domain of F)
%   with parameters A, B.
%
%   If F is an array-valued BNDFUN, then a matrix of coefficients is returned so
%   that F(:,k) = C(N+1,k)*P^{(A,B})_N + ... + C(1,k)*P^{(A,B})_0.
%
% See also CHEBCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

c = jaccoeffs(f.onefun, varargin{:});

end

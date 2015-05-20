function c = jaccoeffs(f, varargin)
%JACCOEFFS   Compute Jacobi series coefficients of a FUN object.
%   C = JACCOEFFS(F, N, A, B) returns the Jacobi series coefficients of a FUN F,
%   so that F = C(N+1)*P^{(A,B)}_N + ... + C(1)*P^{(A,B)}_0, where P^{(A,B)}_k
%   is the degree k Jacobi polynomial (scaled to the domain of F) with
%   parameters A, B.
%
%   If F is an array-valued FUN, then a matrix of coefficients is returned so
%   that F(:,k) = C(N+1,k)*P^{(A,B})_N + ... + C(1,k)*P^{(A,B})_0.
%
% See also CHEBCOEFFS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

c = jaccoeffs(f.onefun, varargin{:});

end

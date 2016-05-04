function c = jaccoeffs(f, n, alp, bet)
%JACCOEFFS   Compute Jacobi series coefficients of a CHEBTECH object.
%   A = JACCOEFFS(F, N, ALPHA, BETA) returns the first N+1 coefficients in the
%   Jacobi series expansion of the CHEBTECH F, so that such that F approximately
%   equals A(1) J_N(x) + ... + A(N) J_1(x) + A(N+1) J_0(x) where J_N(x) denotes
%   the N-th Jacobi polynomial with parameters ALPHA and BETA. A is a column
%   vector. If length(F) < N then the additional entries of A are padded with
%   zeros.
%
%   A = JACCOEFFS(F, ALPHA, BETA) will assume that N = length(F).
%   
%   If F is an array-valued CHEBTECH, then a matrix of coefficients is returned
%   so that F(:,k) = A(1,k)*J_N + ... + A(N,k)*J_0.
%
% See also CHEBCOEFFS, LEGCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 3 )
    % Choose n = length(f)
    bet = alp;
    alp = n;
    n = length(f);
end

% Convert from Chebyshev to Jacobi coefficients.
c = cheb2jac(f.coeffs, alp, bet);

% Truncate / pad to length n:
s = size(c);
if ( s(1) > n )
    c = c(1:n, :);
else
    c = [c ; zeros(n-s(1), s(2))];
end

end

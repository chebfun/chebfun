function D = diffmat( n, m )
%DIFFMAT   Differentiation matrices for ultraspherical spectral method.
%
% D = DIFFMAT( N, M ) returns the differentiation matrix that takes N Chebyshev
% coefficients and returns N C^{(M)} coefficients that represent the derivative
% of the Chebyshev series. Here, C^{(K)} is the ultraspherical polynomial basis with
% parameter K.
%
% This command is a private method.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

%DIFFMAT(N, K, N), computes the kth order US derivative matrix
if ( m > 0 )
    D = spdiags((0:n-1)', 1, n, n);
    for s = 1:m-1
        D = spdiags(2*s*ones(n, 1), 1, n, n) * D;
    end
else
    D = speye(n);
end
end

function D = diffmat(n, m)
%DIFFMAT   Differentiation matrices for ultraspherical spectral method.
%   D = DIFFMAT(N, M) returns the differentiation matrix that takes N Chebyshev
%   coefficients and returns N C^{(M)} coefficients that represent the derivative
%   of the Chebyshev series. Here, C^{(K)} is the ultraspherical polynomial basis
%   with parameter K.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    m = 1;
end

% Create the differentation matrix.
if ( m > 0 )
    D = spdiags((0 : n - 1)', 1, n, n);
    for s = 1:m-1
        D = spdiags(2*s*ones(n, 1), 1, n, n) * D;
    end
else
    D = speye(n);
end

end

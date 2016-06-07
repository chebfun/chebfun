function S = convert(A, K1, K2)
%CONVERT   Conversion (transformation) operator for Ultraspherical method.
%   S = CONVERT(A, K1, K2) returns the operators that maps a vector of C^{(K1)}
%   coefficients to C^{(K2)} coefficients, where C^{(K)} is the ultraspherical
%   polynomial basis with parameter K.
%
%   S = CONVERT(A, 0, K2) returns the operator that maps Chebyshev T
%   coefficients of ultraspherical C^{(K2)} coefficients.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    K2 = A.outputSpace;
end

% Obtaining some useful information:
d = A.domain;
n = A.dimension;
numIntervals = length(d) - 1;

% Find the diagonal blocks:
blocks = cell(numIntervals);
for k = 1:numIntervals
    blocks{k} = ultraS.convertmat(n(k), K1, K2);
end

% Assemble:
S = blkdiag(blocks{:});

end

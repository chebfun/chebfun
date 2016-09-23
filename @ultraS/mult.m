function M = mult(A, f, lambda)
%MULT   Multiplication operator for the ultraspherical spectral method. 
%   M = MULT(A, F, lambda) returns the multiplication operator that represents 
%   u(x) -> F(x)u(x), in the C^{(lambda)} ultraspherical polynomial basis. 
% 
%   If lambda = 0, then the operator is Toeplitz-plus-Hankel-plus-rank-1 and
%   represents multiplication in Chebyshev T coefficients.
%
%   If lambda = 1, then the operator is Toeplitz-plus-Hankel and represents
%   multiplication in Chebyshev U or C^{(1)} coefficients. 
% 
%   If lambda > 1, then the operator does not have any Toeplitz/Hankel structure
%   and is constructed using a three-term recurrence.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Obtaining some useful information:
n = A.dimension;
d = A.domain;
f = restrict(f, d);
numIntervals = length(d) - 1;

% Find the diagonal blocks;
blocks = cell(numIntervals);
for k = 1:numIntervals
    blocks{k} = ultraS.multmat(n(k), f.funs{k}, lambda);
end

% Assemble:
M = blkdiag(blocks{:});

end

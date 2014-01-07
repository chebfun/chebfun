function S = convert( A, K1, K2 )
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
%CONVERT(A, K1, K2), convert C^(K1) to C^(K2)
d = A.domain;
n = A.dimension;
numIntervals = length(d) - 1;
% Find the diagonal blocks.
blocks = cell(numIntervals);
for k = 1:numIntervals
    blocks{k} = convertmat(n(k), K1, K2);
end
% Assemble.
S = blkdiag(blocks{:});
end

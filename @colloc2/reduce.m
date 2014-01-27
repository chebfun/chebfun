function [PA,P] = reduce(disc,blocks)
%REDUCE Dimension reduction for operator matrix. 
%
%   Each block row of the operator DISC.source has an associated dimension
%   reduction to make room for constraints. Given discretized results in BLOCKS,
%   the output PA has one cell per block row, with the resulting projected
%   matrix. The output P has one cell per block row with the projection operator
%   for the row.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

r = sizeReduction(disc.source);
dim = disc.dimension;
PA = cell(size(blocks,1),1);
P = cell(size(blocks,1),1);

for i = 1:size(blocks,1)      % for each block row
    M = cat(2,blocks{i,:});   % combine the blocks for the row into a matrix
    [PA{i},P{i}] = reduceOne(disc,M,dim-r(i),dim);  % do reduction
end

end


function [PA,P] = reduceOne(disc,A,m,n)

% Does reduction for one block row.

if ( m == n )   % do nothing
    PA = A;
    P = eye(size(A,1));
    return
end
        
% Step by intervals in the domain.
domain = disc.domain;
numint = disc.numIntervals;
P = cell(1,numint);

for k = 1:numint
    disc.domain = domain(k:k+1);
    disc.dimension = m(k);
    xOut = points(disc,1);  % projection into 1st kind points
    disc.dimension = n(k);
    xIn = points(disc,2);
    P{k} = barymat(xOut,xIn);
end

P = blkdiag(P{:});
PA = P*A;

end

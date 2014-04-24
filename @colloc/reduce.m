function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Dimension reduction for operator matrix. 
%   Each block row of the operator DISC.source has an associated dimension
%   reduction to make room for constraints. Given discretized results in BLOCKS,
%   the output PA has one cell per block row, with the resulting projected
%   matrix. The output P has one cell per block row with the projection operator
%   for the row.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Setup:
r = sizeReduction(disc.source);
dim = disc.dimension;
dimAdjust = disc.dimAdjust(1,:);
if ( numel(dimAdjust) == 1 )
    dimAdjust = repmat(dimAdjust, 1, size(A, 2));
end
PA = cell(1, size(A, 2));
P = cell(1, size(A, 2));

% Do reduction for each block column:
for i = 1:size(A, 2) 
    [PA{i}, P{i}] = reduceOne(disc, A(:,i), r(i), dim + dimAdjust(i));  
end

% Convert cell arrays to matrices:
P = blkdiag(P{:});
PA = cell2mat(PA);
PS = P;

end


function [PA, P] = reduceOne(disc, A, m, n)
% Does reduction for one block row.

% Step by intervals in the domain.
domain = disc.domain;
numInt = disc.numIntervals;
P = cell(1, numInt);

% Loop through intervals
for k = 1:numInt
    disc.domain = domain(k:(k+1));
    disc.dimension = n(k)-m;
    xOut = equationPoints(disc);
    disc.dimension = n(k);
    [xIn, ~, baryWt] = functionPoints(disc);
    % Store the kth projection matrix in the cell P
    P{k} = barymat(xOut, xIn, baryWt);
end
% Convert the projection matrices P into a blockdiagonal matrix.
P = blkdiag(P{:});

% Project each of the entries of A:
PA = cell(size(A));
for j = 1:numel(A)
    if ( size(P, 2) == size(A{j}, 1) )
        PA{j} = P*A{j};
    else
        PA{j} = A{j};
    end
end
PA = cell2mat(PA);

if ( (m == 0) && (size(A{1}, 2) < sum(n)) )
    % We don't want to project scalars.
    P = eye(size(A, 2));
end

end

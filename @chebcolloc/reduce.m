function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Dimension reduction for operator matrix. 
%   PA = REDUCE(DISC, A) reduces the row dimension of each block column in the
%   cell array A (which is typically a discretization of DISC.SOURCE) so that
%   the reduced discretization, PA, can be formed as a matrix. In particular, PA
%   will have sum(DISC.dimension) rows and sum(cellfun(@(a) size(A, 2), A(1,:))
%   columns.
%
%   [PA, P] = REDUCE(DISC, A) returns also the block-diagonal reduction matrix
%   P. For COLLOC discretizations, this blocks are BARYMAT projections.
%
%   [PA, P, PS] = REDUCE(DISC, A, S) is required for consistency with other
%   opDiscretization reductions. Here S is ignored and PS = P.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Setup:
r = disc.projOrder;
dim = disc.dimension;
dimAdjust = disc.dimAdjust(1,:);
if ( numel(dimAdjust) == 1 )
    dimAdjust = repmat(dimAdjust, 1, size(A, 2));
end
PA = cell(1, size(A, 2));
P = cell(1, size(A, 2));

% Do reduction for each block column:
for k = 1:size(A, 2) 
    [PA{k}, P{k}] = reduceOne(disc, A(:,k), r(k), dim + dimAdjust(k));  
%     [PS{k}, ignored] = reduceOne(disc, S(:,k), r(k), dim + dimAdjust(k));
end

% Convert cell arrays to matrices:
P = blkdiag(P{:});
PA = cell2mat(PA);
% PS = cell2mat(PS);
PS = P;

end

function [PA, P] = reduceOne(disc, A, m, n)
%REDUCEONE   Reduce one block column.
%   [PA, P] = REDUCEONE(DISC, A, M, N) reduces entries of the column cell arrays
%   A from a sum(N)xsum(N) discretization to sum(N-M)xsum(N) version, PA, using
%   the block-projection operator P.

% Step by intervals in the domain.
domain = disc.domain;
numInt = disc.numIntervals;
P = cell(1, numInt);

% Loop through intervals
for k = 1:numInt
    disc.domain = domain(k:(k+1));
    disc.dimension = n(k)-m;
    [xOut, ignored, ignored, tOut] = equationPoints(disc);
    disc.dimension = n(k);
    [xIn, ignored, baryWt, tIn] = functionPoints(disc);
    % Store the kth projection matrix in the cell P
%     P{k} = barymat(xOut, xIn, baryWt);
    P{k} = barymat(xOut, xIn, baryWt, tOut, tIn, 1);
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

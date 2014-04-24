function [PA, P, PS] = reduce(disc, A, S)
%REDUCE Row dimension reduction for operator's matrix. 

%   Each block row of the operator DISC.source has an associated dimension
%   reduction to make room for constraints. Given discretized results in BLOCKS,
%   the output PA has one cell per block row, with the resulting projected
%   matrix. The output P has one cell per block row with the projection operator
%   for the row.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Setup
r = sizeReduction(disc.source);
dim = disc.dimension;
dimAdjust = disc.inputDimensionAdjustment(1,:);

if ( numel(dimAdjust) == 1 )
    dimAdjust = repmat(dimAdjust, 1, size(A, 2));
end

% Outputs will be cells for convenience:
PA = cell(1, size(A, 2));
P = cell(1, size(A, 2));

for i = 1:size(A, 2) % Do reduction for each block column:
    [PA{i}, P{i}] = reduceOne(disc, A(:,i), r(i), dim + dimAdjust(i));  
end

PS = S;
for j = 1:size(S, 1)
    for k = 1:size(S, 2)
        if ( size(A{j,k}, 1) > 1 )
            PS{j,k} = P{k}*S{j,k};
        end
    end
end

% Convert cell arrays to matrices:
P = blkdiag(P{:});
PA = cell2mat(PA);
PS = cell2mat(PS);

end


function [PA, P] = reduceOne(disc, A, m, n)
% TODO: What do the input variables stand for? DISC is clear
% from above, presumably, A is a block-row. What does m do?
% AB, 1/3/14.
dom = disc.domain;

% chop off some rows and columns
% TODO: It's not really obvious what's going on here. AB, 1/3/14.
v = [];
nn = cumsum([0 n]);
P = {};
for k = 1:numel(dom) - 1
    v = [v nn(k) + n(k) - (0:m-1)];
end
P = speye(sum(n));
P(v,:) = [];

% Convert the projection matrices P into a blockdiagonal matrix.
PA = A;
for j = 1:numel(PA)
    if ( size(P, 2) == size(A{j}, 1) )
        PA{j}(v.', :) = [];
    else
        PA{j} = A{j};
    end
end
PA = cell2mat(PA);

if ( m == 0 && size(A{1},2) < sum(n) )
    % We don't want to project scalars.
    P = speye(size(A, 2));
end

end

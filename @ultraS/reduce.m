function [PA, P] = reduce(disc, blocks)
%REDUCE Row dimension reduction for operator's matrix. 

%   Each block row of the operator DISC.source has an associated dimension
%   reduction to make room for constraints. Given discretized results in BLOCKS,
%   the output PA has one cell per block row, with the resulting projected
%   matrix. The output P has one cell per block row with the projection operator
%   for the row.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

r = sizeReduction(disc.source);
dim = disc.dimension;
space = disc.inputDimension(1,:);

PA = cell(1, size(blocks, 1));
P = cell(1, size(blocks, 1));
% Loop through the block-rows of the operator.

for i = 1:size(blocks, 2)       % for each block column
    [PA{i}, P{i}] = reduceOne(disc, blocks(:,i), r(i), dim+space(i));  % do reduction
end

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
    v = [v nn(k) + n(k) - (1:m)];
end
P = eye(sum(n));
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
    P = eye(size(A, 2));
end

end

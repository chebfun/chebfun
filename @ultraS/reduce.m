function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Row dimension reduction for operator's matrix. 
%   Each block row of the operator DISC.source has an associated dimension
%   reduction to make room for constraints. Given discretized results in BLOCKS,
%   the output PA has one cell per block row, with the resulting projected
%   matrix. The output P has one cell per block row with the projection operator
%   for the row.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Setup
r = sizeReduction(disc.source);
dim = disc.dimension;
dimAdjust = disc.dimAdjust(1,:);
if ( numel(dimAdjust) == 1 )
    dimAdjust = repmat(dimAdjust, 1, size(A, 2));
end
PA = cell(1, size(A, 2));
P = cell(1, size(A, 2));
PS = cell(1, size(A, 2));

for i = 1:size(A, 2) % Do reduction for each block column:
    [PA{i}, P{i}, PS{i}] = reduceOne(disc, A(:,i), S(:,i), r(i), dim + dimAdjust(i));  
end

% Convert cell arrays to matrices:
P = blkdiag(P{:});
PA = cell2mat(PA);
PS = cell2mat(PS);

end


function [PA, P, PS] = reduceOne(disc, A, S, m, n)
% TODO: What do the input variables stand for? DISC is clear
% from above, presumably, A is a block-row. What does m do?
% AB, 1/3/14.
dom = disc.domain;

% Projection matrix for US removes the last m coeffs:
P = speye(sum(n));
v = []; % Row indicies which are to be removed by projection.
nn = cumsum([0 n]);
for k = 1:numel(dom) - 1
    v = [v nn(k) + n(k) - (0:m-1)];
end
P(v,:) = [];

% Project each component of A and S:
PA = A; 
PS = S;
for j = 1:numel(PA)
    if ( size(P, 2) == size(A{j}, 1) )
        PA{j}(v.', :) = [];
        PS{j}(v.', :) = [];
    else
        PA{j} = A{j};
        PS{j} = S{j};
    end
end
PA = cell2mat(PA);
PS = cell2mat(PS);

if ( (m == 0) && (size(A{1}, 2) < sum(n)) )
    % We don't want to project scalars.
    P = eye(size(A, 2));
end

end

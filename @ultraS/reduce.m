function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Dimension reduction for operator matrix. 
%   PA = REDUCE(DISC, A) reduces the row dimension of each block column in the
%   cell array A (which is typically a discretization of DISC.SOURCE) so that
%   the reduced discretization, PA, can be formed as a matrix. In particular, PA
%   will have sum(DISC.dimension) rows and sum(cellfun(@(a) size(A, 2), A(1,:))
%   columns.
%
%   [PA, P] = REDUCE(DISC, A) returns also the block-diagonal reduction matrix
%   P. For ultraS discretizations, this is a sparse identity matrix with some
%   missing rows.
%
%   [PA, P, PS] = REDUCE(DISC, A, S) returns also a matrix of the projected
%   conversion operator cell array S.

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
%REDUCEONE   Reduce one block column.
%   [PA, P, PS] = REDUCEONE(DISC, A, S, M, N) reduces entries of the column cell
%   arrays A and S from sum(N)xsum(N) discretizations to sum(N-M)xsum(N)
%   versions (PA and PS, respectively) using the block-projection operator P.

% Projection matrix for US removes the last m coeffs:
P = speye(sum(n));
v = []; % Row indicies which are to be removed by projection.
nn = cumsum([0 n]);
for k = 1:disc.numIntervals
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

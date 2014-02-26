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
PA = cell(size(blocks,1),1);
P = cell(size(blocks,1),1);
for i = 1:size(blocks,1)
    M = cat(2, blocks{i,:});
    [PA{i}, P{i}] = reduceOne(disc, M, dim - r(i));
end

end


function [A, P] = reduceOne(disc, A, m)

dom = disc.domain;
n = disc.dimension;
% chop off some rows and columns
v = [];
nn = cumsum([0 n]);
P = eye(size(A,1));
for k = 1:numel(dom)-1
    v = [v m(k) + nn(k) + (1:(n(k)-m(k)))];
end
A(v.', :) = [];
P(v.', :) = [];
end

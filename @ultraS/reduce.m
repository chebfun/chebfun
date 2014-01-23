function [B, P] = reduce(disc, blocks)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

r = sizeReduction(disc.source);
dim = disc.dimension;
B = cell(size(blocks,1),1);
P = cell(size(blocks,1),1);
for i = 1:size(blocks,1)
    M = cat(2, blocks{i,:});
    [B{i}, P{i}] = reduceOne(disc, M, dim - r(i));
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

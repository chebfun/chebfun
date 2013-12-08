function [B, P] = reproject(disc, blocks)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
reduce = disc.source.sizeReduction;
dim = disc.dimension;
B = cell(size(blocks,1),1);
P = cell(size(blocks,1),1);
for i = 1:size(blocks,1)
    M = cat(2, blocks{i,:});
    [B{i}, P{i}] = resize(disc, M, dim - reduce(i));
end
end

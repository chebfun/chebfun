function [A, P] = resize(disc, A, m)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
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

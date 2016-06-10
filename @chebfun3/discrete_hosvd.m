function [core, U1, U2, U3] = discrete_hosvd(T)
% Higher order SVD of a discrete tensor of order 3

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Mode-1 unfolding of T:
T1 = chebfun3.unfold(T, 1); 
[U1,~,~] = svd(T1, 'econ');

% Mode-2 unfolding of T:
T2 = chebfun3.unfold(T, 2);
[U2,~,~] = svd(T2, 'econ');

% Mode-3 unfolding of T:
T3 = chebfun3.unfold(T, 3);
[U3,~,~] = svd(T3, 'econ');

core = chebfun3.txm(chebfun3.txm(chebfun3.txm(T, U1', 1), U2', 2), U3', 3);

end
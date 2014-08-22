function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Dimension reduction for operator matrix. 
%   Does not do anything for FOURCOLLOC objects.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

PA = cell2mat(A);
P = eye(length(PA));
PS = [];

end

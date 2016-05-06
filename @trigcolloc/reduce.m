function [PA, P, PS] = reduce(disc, A, S)
%REDUCE   Dimension reduction for operator matrix.
%   PA = REDUCE(DISC, A) convert the cell array A (which is typically a 
%   discretization of DISC.SOURCE) to a matrix.
%
%   [PA, P, PS] = REDUCE(DISC, A, S) is required for consistency with other
%   OPDISCRETIZATION reductions. Here S is ignored, P is the identity and
%   PS = P.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert A to a matrix.
PA = cell2mat(A);

% No projection: the matrix P is the identity.
P = speye(length(PA));
PS = P;

end

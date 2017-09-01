function u = partition(disc, values)
%OPDISCRETIZATION.PARTITION   Partition values to appropriate variables.
%   U = OPDISCRETIZATION.PARTITION(DISC, VALUES) will, given a vector or
%   matrix (columnwise) VALUES of values corresponding to all the discretized
%   variables and scalars in a system DISC, convert to a cell-valued partition
%   of individual variables in the system. I.e., deduce the variable boundaries
%   within the discretization.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Which variables are functions (as opposed to scalars)?
isFun = isFunVariable(disc.source);

% Allocate the discretization size to each function variable.
m = ones(size(isFun));
m(isFun) = sum(disc.dimension);

% Do the partition.
u = mat2cell(values, m, size(values, 2));

end

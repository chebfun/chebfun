function u = partition(disc,values)

% Given a vector or matrix (columnwise) of values corresponding to all the
% discretized variables and scalars, convert to a cell-valued partition of
% individual variables. 

% Which variables are functions (as opposed to scalars)?
isFun = isFunVariable(disc.source);

% Allocate the discretization size to each function variable.
m = ones(size(isFun));
m(isFun) = sum(disc.dimension);

% Do the partition.
u = mat2cell(values, m, size(values,2));

end

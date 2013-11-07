function [A, b, dom] = linSystem(L, f, disc)

% Domain needs to have the union of all breakpoints:
Lblocks = L.operator.blocks;
fblocks = f.blocks;
domain = chebmatrix.mergeDomains( [Lblocks, fblocks] );
disc.domain = domain;

% Discretize the operator and the RHS:
A = matrix(disc);
b = rhs(disc,f);

end

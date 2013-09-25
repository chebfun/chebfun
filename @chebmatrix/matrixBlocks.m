function A = matrixBlocks(L,dim,varargin)
% matrixBlocks(A,DIM) returns a cell array in which each block is represented by its
% DIM-dimensional discretization.
%
% matrixBlocks(A,DIM,TYPE) uses the discretization type TYPE.

p = inputParser;
addOptional(p,'domain',domain(L),@isnumeric);
addOptional(p,'matrixType',linop.defaultDiscretization,@(x) isa(x,'function_handle'));
parse(p,varargin{:});


data = L.blocks;
[m,n] = size(L);
dom = p.Results.domain;

% Discretize each block in place.
A = cell(m,n);
for i = 1:m
    for j = 1:n
        item = data{i,j};
        if isa(item,'linop')
            A{i,j} = matrix(item,dim,dom,p.Results.matrixType);
        elseif isa(item,'chebfun')
            x = chebpts(dim, dom);
            A{i,j} = item(x);
        else   % scalar
            A{i,j} = item;
        end
    end
end
end

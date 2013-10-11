function A = matrixBlocks(L,dim,varargin)
% matrixBlocks(A,DIM) returns a cell array in which each block is represented by its
% DIM-dimensional discretization.
%
% matrixBlocks(A,DIM,TYPE) uses the discretization type TYPE.

p = inputParser;
addOptional(p,'domain',L.domain,@isnumeric);
addOptional(p,'matrixType',linBlock.defaultDiscretization,@(x) isa(x,'function_handle'));
parse(p,varargin{:});

data = L.blocks;
[m,n] = size(L);
dom = p.Results.domain;

% Discretize each block in place.
A = cell(m,n);
dummy = p.Results.matrixType([]);
for i = 1:m
    for j = 1:n
        item = data{i,j};
        if isa(item,'linBlock')
            A{i,j} = discretize(item,dim,dom,p.Results.matrixType);
        elseif isa(item,'chebfun')
            A{i,j} = dummy.discretizeFunction(item,dim,dom);
        else   % scalar
            A{i,j} = item;
        end
    end
end
end

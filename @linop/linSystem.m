function [A,b,dom] = linSystem(L,f,dim,matrixType)

if nargin < 4
    matrixType = linBlock.defaultDiscretization;
end

if isa(f,'chebfun')
    f = chebmatrix({f});
end

% Domain needs to have the union of all breakpoints.
Lblocks = L.operator.blocks;
fblocks = f.blocks;
dom = chebmatrix.mergeDomains( {Lblocks{:}, fblocks{:}} );
L.operator.domain = dom;
f.domain = dom;

% Update domain of constraints.
L.constraint.operator.domain = dom;

L = appendContinuity(L);

% This is needed when generating chebpts for discretization of a chebfun.
if (length(dim)==1)
    dim = repmat(dim,1,length(dom)-1);
end

Ablocks = discretizeBlocks(L.operator,dim,dom,matrixType);
bblocks = discretizeBlocks(f,dim,dom,matrixType);
% Resize the blocks and append the constraints.
[m,n] = size(L);
rows = cell(m+1,1);

% Resize the operator rows according to differential order.
dummy = matrixType([]);
d = getDownsampling(L);
for i = 1:m
    M = cat(2,Ablocks{i,:},bblocks{i});
    if ( ~isnan(d(i)) ) && ( d(i) > 0 )
        M = dummy.resize( M, dim-d(i), dim, dom );
    end
    rows{i} = M;
end

% Append the discrete constraints.
bcOp = L.constraint.operator;
bcVal = L.constraint.values;
rows{m+1} = [ discretize(bcOp,dim,dom,matrixType), bcVal ];
aug = cell2mat(rows);
A = aug(:,1:end-1);
b = full(aug(:,end));

end

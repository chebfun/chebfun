function [A,b,dom] = linSystem(L,f,dim,matrixType)

if nargin < 4
    matrixType = linop.defaultDiscretization;
end

if isa(f,'chebfun')
    f = chebmatrix({f});
end

% Domain needs to have the union of all breakpoints.
dom = chebmatrix.mergeDomains( {L.blocks{:}, f.blocks{:}} );
L.domain = dom;
f.domain = dom;

% Update domains of constraints:
for k = 1:numel(L.constraints)
    L.constraints(k).op.domain = dom;
end

L = appendContinuity(L);

% This is needed when generating chebpts for discretization of a chebfun.
if (length(dim)==1)
    dim = repmat(dim,1,length(dom)-1);
end

Ablocks = matrixBlocks(L,dim,dom,matrixType);
bblocks = matrixBlocks(f,dim,dom,matrixType);

% Resize the blocks and append the constraints.
[m,n] = size(L);
rows = cell(m+numbc(L),1);

% Resize the operator rows according to differential order.
dummy = matrixType([]);
d = getDownsampling(L);
for i = find( ~isnan(d) )
    M = cat(2,Ablocks{i,:},bblocks{i});
    rows{i} = dummy.resize( M, dim-d(i), dim, dom );
end

% Append the discrete constraints.
bc = L.constraints;
for i = 1:length(bc)
    M = [ matrix(bc(i).op,dim,dom,matrixType), bc(i).value ];
    rows{m+i} = M;
end

aug = cell2mat(rows);
A = aug(:,1:end-1);
b = aug(:,end);

end

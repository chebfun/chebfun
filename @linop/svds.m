function [U, S, V] = svds(L, k, bcType)
%SVDS  Find some singular values and vectors of a compact LINOP.
%   SVDS of a LINOP is currently not supported.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% construct adjoint
Lstar = adjoint(L,bcType);

% initialize superL to 0
[m n] = size(L);
dom = L.domain;
nm = n + m;
superL = linop(mat2cell(zeros(nm),ones(1,nm),ones(1,nm)));
z = chebfun(0,dom);
for ii = 1:nm
    for jj = 1:nm
        superL.blocks{ii,jj} = operatorBlock.mult(z,z.domain);
    end
end
superL.domain = L.domain;

% fill the appropriate blocks
for ii = 1:m
    for jj = 1:n
        superL.blocks{n+ii,jj} = L.blocks{ii,jj};
    end
end
for ii = 1:n
    for jj = 1:m
        superL.blocks{ii,n+jj} = Lstar.blocks{ii,jj};
    end
end

% get constraints
C = L.constraint;
Cstar = Lstar.constraint;

% get funs
Cfuns = C.functional; nc = size(Cfuns,1);
Csfuns = Cstar.functional; ncs = size(Csfuns,1);

% make superC
dor = max(max(L.diffOrder));
superC.functional = [];
z = functionalBlock.zero(dom);
for ii = 1:nc
    row = [];
    for jj = 1:n
        row = [row,Cfuns{ii,jj}];
    end
    for jj = 1:m
        row = [row,z];
    end
    superC.functional = [superC.functional;row];
end
for ii = 1:ncs
    row = [];
    for jj = 1:n
        row = [row,z];
    end
    for jj = 1:m
        row = [row,Csfuns{ii,jj}];
    end
    superC.functional = [superC.functional;row];
end
superC.values = zeros(dor*nm,1);

% finish superL
superL.constraint = superC;

% count number of null vectors for L and Lstar
nulL = dor*n-nc

% set number of singular values
nsvals = 2*k-abs(nulL) 

% call linop/eigs
pref = cheboppref();
if ( strcmp(bcType,'periodic') )
    pref.discretization = @trigcolloc;
else
    pref.discretization = @chebcolloc2;
end
[ Q, D ] = eigs( superL, nsvals, [], pref, 'rayleigh' );

% make sure singular values are real
if ( any(imag(diag(D)) ~= 0) )
    error('CHEBFUN:LINOP:svds:real', ...
        'Computed singular values are not strictly real.');
end

% sort
[D,id] = sort(diag(D),'descend');
Q = Q(:,id);

% trim 
S = diag(D(1:k));
Q = Q(:,1:k);

% rescale singular vectors
V = Q(1:n,:); nrmV = sqrt(diag(V'*V)); 
nrmV( nrmV < pref.bvpTol ) = 1; V = V*diag(1./nrmV);
U = Q(n+1:end,:); nrmU = sqrt(diag(U'*U));
nrmU( nrmU < pref.bvpTol ) = 1; nrmU = 1./nrmU;
U = U*diag(nrmU);

end

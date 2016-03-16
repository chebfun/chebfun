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

% call linop/eigs
pref = cheboppref();
if ( strcmp(bcType,'periodic') )
    pref.discretization = @trigcolloc;
else
    pref.discretization = @chebcolloc2;
end
[ Q, D ] = eigs( superL, 2*k, [], pref );

% make sure singular values are real
if ( any(imag(diag(D)) ~= 0) )
    error('CHEBFUN:LINOP:svds:real', ...
        'Computed singular values are not strictly real.');
end

% sort
[D,id] = sort(diag(D),'descend');
Q = Q(:,id);

% one step of RQI
I = operatorBlock.eye(superL.domain);
for ii = 1:2*k
    K = superL;
    for jj = 1:size(K,1)
        K.blocks{jj,jj} = K.blocks{jj,jj} - D(ii)*I;
    end
    Q(:,ii) = linsolve(K,Q(:,ii),pref);
    Q(:,ii) = Q(:,ii)/norm(Q(:,ii));
end
D = diag(Q'*(superL*Q))

% check for zero singular values
% a singular value is considered 0 if it is less than pref.bvpTol
nulls = abs(D) < pref.bvpTol;
ids = (1:length(nulls))';
zid = max(ids(nulls)); 
nzs = sum([nulls;0]); 
if ( ~isempty(zid) )
    ids = ids(zid-k+1:zid);
else
    ids = ids(1:k);
end

% trim 
S = diag(D(ids));
Q = Q(:,ids);

% rescale singular vectors
V = Q(1:n,:); nrmV = sqrt(diag(V'*V)); 
nrmV( nrmV == 0 ) = 1; V = V*diag(1./nrmV);
U = Q(n+1:end,:); nrmU = sqrt(diag(U'*U));
nrmU( nrmU == 0 ) = 1; nrmU = 1./nrmU;
if ( nzs > 0 )
    nrmU(end-nzs+1:end) = 0;
end
U = U*diag(nrmU);

% try to simply null vectors
nulV = V(:,end-nzs+1:end);
nulV = squish(nulV,pref);
V(:,end-nzs+1:end) = nulV;

end

function Vnew = squish(V,pref)
% this function attempts construct a degree-graded 
% set of vectors Vnew from the input vectors V

% convert V to quasimatrix
V = quasi2cheb(quasimatrix(V));
V = simplify(V,1e-8);

% extract coefficients
Vcfs = chebcoeffs(V);

% reduce using partial pivoting
Vcfs = rref(fliplr(Vcfs'));
Vcfs = fliplr(Vcfs)';
V = chebfun(Vcfs,domain(V),'coeffs');

% reorthogonalize
[V,~] = qr(V);

% set Vnew
Vnew = chebmatrix(V);

end

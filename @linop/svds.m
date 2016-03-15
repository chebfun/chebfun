function [U, S, V] = svds(L, k, bcType)
%SVDS  Find some singular values and vectors of a compact LINOP.
%   SVDS of a LINOP is currently not supported.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% construct adjoint
Lstar = adjoint(L,bcType);

% initialize superL to 0
[m n] = size(L)
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

% make superC
superC.functional = cell(nm,nm)
z = functionalBlock.zero(dom);
for ii = 1:m
    for jj = 1:n
        superC.functional{ii,jj} = C.functional{ii,jj};
    end
    for jj = 1:m
        superC.functional{ii,n+jj} = z;
    end
end
for ii = 1:n
    for jj = 1:n
        superC.functional{m+ii,jj} = z;
    end
    for jj = 1:m
        superC.functional{m+ii,n+jj} = Cstar.functional{ii,jj};
    end
end
superC.functional
superC.values = zeros(nm,1);

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

% sort and trim
[D,id] = sort(diag(D),'descend');
Q = Q(:,id);
S = diag(D(1:k));
Q = Q(:,1:k);

% rescale singular vectors
V = Q(1:n,:); nrmV = sqrt(diag(V'*V)); 
nrmV( nrmV == 0 ) = 1; V = V*diag(1./nrmV);
U = Q(n+1:end,:); nrmU = sqrt(diag(U'*U));
nrmU( nrmU == 0 ) = 1; U = U*diag(1./nrmU);

end

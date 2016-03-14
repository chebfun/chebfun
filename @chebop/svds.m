function [U, S, V] = svds(L, k, sigma)
%SVDS  Find some singular values and vectors of a compact linear CHEBOP.
%   SVDS of a CHEBOP is currently not supported.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%error('CHEBFUN:CHEBOP:svds:nosupport', ...
%    'CHEBOP/SVDS() is not currently supported.');

% default outputs
U = [];
S = [];
V = [];

% Linearize
K = linearize(L);

% construct adjoint
Kstar = adjoint(K,getBCType(L));

% initialize superK to 0
[m n] = size(K);
nm = n + m;
superK = linop(mat2cell(zeros(nm),ones(1,nm),ones(1,nm)));
z = chebfun(0,K.domain);
for ii = 1:nm
    for jj = 1:nm
        superK.blocks{ii,jj} = operatorBlock.mult(z,z.domain);
    end
end
superK.domain = K.domain;

% fill the appropriate blocks
for ii = 1:m
    for jj = 1:n
        superK.blocks{n+ii,jj} = K.blocks{ii,jj};
    end
end
for ii = 1:n
    for jj = 1:m
        superK.blocks{ii,n+jj} = Kstar.blocks{ii,jj};
    end
end

% get constraints
C = K.constraint;
Cstar = Kstar.constraint;

% make superC
superC.functional = [ C.functional, 0*Cstar.functional; ... 
                      0*C.functional, Cstar.functional ];
superC.values = 0*[ C.values; Cstar.values ];

% finish superK
superK.constraint = superC;

% call linop/eigs
pref = cheboppref();
pref.discretization = @chebcolloc2;
[ Q, D ] = eigs( superK, 2*k, sigma, pref );

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

function [Q, R] = qr(A,econ)
%QR   QR factorization of an array-valued quasimatrix.
%   [Q,R] = QR(A) or QR(A, 0), where A is a column CHEBFUN with n columns,
%   produces a column CHEBFUN Q with n orthonormal columns and an n x n upper
%   triangular matrix R such that A = Q*R.
%
%   This algorithm used is described in L.N. Trefethen, "Householder
%   triangularization of a quasimatrix", IMA J. Numer. Anal. (30), 887-897
%   (2010).
%
% See also SVD, MRDIVIDE, RANK.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check inputs
if ( ( nargin > 2 ) || ( nargin == 2 && econ ~= 0 ) )
    error('CHEBFUN:qr:twoargs',...
        'Use qr(A) or qr(A,0) for QR decomposition of an array-valued CHEBFUN.');
end
if ( A.isTransposed )
    error('CHEBFUN:qr:transpose',...
        'CHEBFUN QR works only for column CHEBFUN objects.')
end
if ( ~all(isfinite(A.domain)) )
    error('CHEBFUN:QR:infdomain', ...
        'CHEBFUN QR does not support unbounded domains.');
end

% No breakpoints = easy case:
if ( numel(A.funs) == 1 )
    [Q, R] = qr(A.funs{1});
    Q = chebfun({Q});
    return
end

error
% TODO: The type of indexing used below does not work for array-valued CHEBFUNs.
% TODO: Imple,ent more efficiently.

% Get some useful values
n = size(A, 2);
R = zeros(n);
tol = epslevel(A)*vscale(A);
a = A.domain(1);
b = A.domain(end);

% Set up target quasimatrix E with orthonormal columns:
E = legpoly(0:n-1, [a, b], 'norm');

% Householder triangularization:
V = chebfun();                        % cols of V will store Househ. vectors
for k = 1:n
    I = 1:k-1; J = k+1:n;             % convenient abbreviations
    e = E(k);                         % target for this reflection
    x = A(k);                         % vector to be mapped to s*r*e
    ex = e'*x; aex = abs(ex);
    if ( aex == 0 )
        s = 1; 
    else
        s = -ex/aex; 
    end
    e = s*e; E(k) = e;                % adjust e by sign factor
    r = norm(x); R(k,k) = r;          % diagonal entry r_kk
    v = r*e - x;                      % vector defining reflection
    if ( k > 1 )
        v = v - E(I)*(E(I)'*v);       % improve orthogonality
    end
    nv = norm(v);
    if ( nv < tol*max(vscale(x), vscale(e)) )
        v = e;
    else
        v = v/nv;
    end
    V(k) = v;                         % store this Householder vector
    if ( k < n )
        A(J) = A(J)-2*v*(v'*A(J));    % apply the reflection to A
        rr = e'*A(J); R(k,J) = rr;    % kth row of R
        A(J) = A(J) - e*rr;           % subtract components in direction e
    end
end

% Form the quasimatrix Q from the Householder vectors:
Q = E;
for k = n:-1:1
    v = V(k);
    J = k:n;
    w = v'*Q(J);
    Q(J) = Q(J) - 2*v*w;
end

end

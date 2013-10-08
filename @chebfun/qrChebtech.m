function [Q, R] = qrChebtech(A, econ)
%QR   QR factorization of an array-valued quasimatrix.
%   [Q, R] = QR(A) or QR(A, 0), where A is a column CHEBFUN with n columns,
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

% [TODO]: This implementation is fast, but relies on the fact that everything is
% a CHEBTECH at heart.

% Check inputs
if ( ( nargin > 2 ) || ( nargin == 2 && econ ~= 0 ) )
    error('CHEBFUN:qr:twoargs',...
      'Use qr(A) or qr(A, 0) for QR decomposition of an array-valued CHEBFUN.');
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

A = simplify(A);

% Get some useful values
numCols = min(size(A));
R = zeros(numCols);
tol = epslevel(A)*vscale(A);
a = A.domain(1);
b = A.domain(end);
kind = 2;
ends = A.domain;
iends = ends(2:end-1).';
numFuns = length(ends)-1;

% Get the sizes of the funs in the columns of A, keeping in mind that we
% will have to multiply with the columns of E and the norm of A's columns
sizes = zeros(numFuns, 1);
for j = 1:numFuns
    sizes(j) = 2*max(length(A.funs(j)), numCols);
    A.funs{j}.onefun = prolong(A.funs{j}.onefun, sizes(j));
end
inds = [ 0 ; cumsum(sizes) ];

% Create the chebyshev nodes and quadrature weights
[pts, w] = chebpts( sizes , ends , kind );

    % ------------------------------------
    % Define the inner product as a nested function:
    function res = innerProduct(f, g)
        res = w * (conj(f) .* g);
    end
    % ------------------------------------

% Make the discrete Analog of A
dA = get(A, 'values');
dA = cat(1, dA{:});

% Generate a discrete E (Legendre-Chebyshev-Vandermonde matrix) directly:
xx = 2*(pts - (a + b)/2) / (b - a); % Unscale the CHEBPTS().
dE = ones(size(dA));
dE(:,2) = xx;
for k = 3:numCols % Recurrence relation:
    dE(:,k) = ( (2*k - 3)*xx.*dE(:,k-1) - (k - 2)*dE(:,k-2) ) / (k - 1);
end
% Scaling:
for k = 1:numCols
    dE(:,k) = dE(:,k) * sqrt( (2*k - 1) / (b - a) );
end

% Pre-allocate the matrix V:
V = zeros(size(dA));


% Do the QR-thing, but with the discretized values:
for k = 1:numCols

    % Indices of the previous and following columns:
    I = 1:k-1; J = k+1:numCols;
    
    % Scale:
    scl = max(norm(dE(:,k),inf), norm(dA(:,k), inf));

    % Multiply the kth column of A with the basis in E:
    ex = innerProduct(dE(:,k), dA(:,k) );
    aex = abs(ex);

    % Adjust the sign of the kth column in E:
    if ( aex < eps*scl )
        s = 1; 
    else
        s = -sign(ex/aex);
    end
    dE(:,k) = dE(:,k) * s;

    % Compute the norm of the kth column of A:
    r = sqrt(innerProduct(dA(:,k), dA(:,k)) );
    R(k,k) = r;

    % Compute the reflection v:
    v = r*dE(:,k) - dA(:,k);
    % Make it more orthogonal:
    for i = I
        ev = innerProduct( dE(:,i) , v );
        v = v - dE(:,i)*ev;
    end
    % Normalize:
    nv = sqrt( innerProduct( v , v ) );
    if ( nv < tol*scl )
       v = dE(:,k); 
    else
       v = v / nv; 
    end
    % Store:
    V(:,k) = v;

    % Subtract v from the remaining columns of A:
    for j = J
        av = innerProduct(v, dA(:,j));
        dA(:,j) = dA(:,j) - 2*v*av;
        rr = innerProduct( dE(:,k) , dA(:,j) );
        dA(:,j) = dA(:,j) - dE(:,k)*rr;
        R(k,j) = rr;
    end

end

% Form a discrete Q from the columns of V:
dQ = dE;
for k = numCols:-1:1
    for j = k:numCols
        vq = innerProduct(V(:,k), dQ(:,j));
        dQ(:,j) = dQ(:,j) - 2*(V(:,k)*vq);
    end
end

% Construct a CHEBFUN from the discrete values:
dQ = mat2cell(dQ, sizes, numCols);
Q = chebfun(dQ, ends);

end
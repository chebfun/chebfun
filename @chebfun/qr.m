function [Q, R] = qr(A, econ, flag)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%  If A contains only a single breakpoint, then FUN/QR is used directly. If A
%  has multiple pieces but each of these are simple CHEBTECH objects, then
%  QRCHEBTECH() is called. This violates OOP principles, but is _much_ more
%  efficient.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [TODO]: The third input argument is for testing qrChebfun. It should be
% removed and a test added which requires this code (e.g., 

% Check inputs
% if ( ( nargin > 2 ) || ( nargin == 2 && econ ~= 0 ) )
if ( ( nargin == 2 && econ ~= 0 ) ) % TODO: Replace above once flag is removed.
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
econ = 0;

% No breakpoints = easy case:
if ( numel(A.funs) == 1 )
    [Q, R] = qr(A.funs{1});
    Q = chebfun({Q});
    return
end

% If all the FUN objects are on BNDFUNs and their techs are CHEBTECHs, we can
% use a much more efficient approach. Note, this completely violates OOP
% principles, but the performance gain is worth it.
isSimple = all(cellfun(@(f) isa(f.onefun, 'chebtech'), A.funs));
if ( isSimple && nargin < 3 )
    [Q, R] = qrChebtech(A, econ);
    return
end

% Use L.N. Trefethen, "Householder  triangularization of a quasimatrix" in the
% continuous setting:
[Q, R] = qrChebfun(A, econ);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R] = qrChebfun(A, varargin)

% Get some useful values
n = size(A, 2);
R = zeros(n);
tol = epslevel(A)*vscale(A);
a = A.domain(1);
b = A.domain(end);

% Simplify A:
A = simplify(A);

% Set up target quasimatrix E with orthonormal columns:
E = legpoly(0:n-1, [a, b], 'norm');
[A, E] = overlap(A, E);

% Householder triangularization:
V = cell(1, n);                       % cols of V will store Househ. vectors
for k = 1:n
    I = 1:k-1; J = k+1:n;             % convenient abbreviations
    e = extractColumns(E, k);         % target for this reflection
    x = extractColumns(A, k);         % vector to be mapped to s*r*e
    
    ex = innerProduct(e, x);
    aex = abs(ex);
    if ( aex == 0 )
        s = 1; 
    else
        s = -ex/aex; 
    end
    e = s*e;                          % adjust e by sign factor
    E = assignColumns(E, k, e);
    
    r = norm(x); R(k,k) = r;          % diagonal entry r_kk
    v = r*e - x;                      % vector defining reflection
    if ( k > 1 )
        Ei = extractColumns(E, I);
        c = innerProduct(Ei, v);
        v = v - Ei*c;                 % improve orthogonality
    end
    nv = norm(v);
    if ( nv < tol*max(vscale(x), vscale(e)) )
        v = e;
    else
        v = v/nv;
    end
    
    if ( k < n )
        Aj = extractColumns(A, J);
        c = 2*innerProduct(v, Aj);
        Aj = Aj - v*c;                % apply the reflection to A
        rr = innerProduct(e, Aj); 
        R(k,J) = rr;                  % kth row of R
        Aj = Aj - e*rr;               % subtract components in direction e
        A = assignColumns(A, J, Aj);  % assign back to A
    end
    
    V{k} = v;                         % store this Householder vector
end

% Form the quasimatrix Q from the Householder vectors:
Q = E;
for k = n:-1:1
    v = V{k};
    J = k:n;
    Qj = extractColumns(Q, J);
    c = 2*(v'*Qj);
    Qj = Qj - v*c;
    Q = assignColumns(Q, J, Qj);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R] = qrChebtech(A, varargin)
% This implementation is fast, but relies on the fact that everything is a
% CHEBTECH on a bounded domainat heart.

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

% Simplify A:
A = simplify(A);

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
    ex = innerProduct(dE(:,k), dA(:,k));
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
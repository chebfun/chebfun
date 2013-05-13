function [f, R, E] = qr2(f, flag)
%QR   QR factorisation of an array-valued CHEBTECH.
%   [Q, R] = QR2(F) returns a QR factorisation of F such that F = Q*R, where the
%   CHEBTECH Q is orthogonal (with respect to the continuous L^2 norm on [-1,1])
%   and of the same size as F and R is an m x m upper-triangular matrix when F
%   has m columns. It is an implementation of the algorithm described in [1].
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.
%
% See also CHEBTECH/QR.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This version of QR is that suggested by Trefethen in [1]. Typically this will
% be slower than CHEBTECH/QR(), but should be more stable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with empty case:
if ( isempty(f) )
    R = [];
    return
end

% Grab the size of f:
m = size(f, 2);

% If f has only one column we simply scale it.
if ( m == 1 )
    R = sqrt(innerProduct(f, f));
    f = f./R;
    E = 1;
    return
end

% Specify a tolerance:
pref = chebtech.pref;
tol = pref.chebtech.eps;

% Simplify so that we don't do any extra work: (QR is O(m*n^2)? :/ )
f = simplify(f, pref);
n = size(f, 1);

% Create the Chebyshev nodes and quadrature weights of double the length:
% (Note: we double the size so that sum(f*f) is exact for f in P_{n}.)
x = f.chebpts(2*n);
w = f.quadwts(2*n);

    % Define the inner product as a nested function:
    function res = innerprod(a, b)
        res = w*(conj(a).*b);
    end

% Make the matrix of function values on a doubled grid so that Q*R=A:
A = get(prolong(f, 2*n), 'values');

% Generate the Legendre-Vandermonde matrix E via the standard recurrence:
E = ones(2*n, m);
E(:,2) = x;
for k = 3:m
    E(:,k) = ((2*k-3)*x.*E(:,k-1) - (k - 2)*E(:,k-2)) / (k - 1);
end
for k = 1:m
    E(:,k) = E(:,k) * sqrt(k - .5);
end

% Pre-allocate memory for R and V:
R = zeros(m);
V = zeros(2*n, m);

% Discretised version of code from Trefethen's paper:
for k = 1:n
    
    % Indices of the previous and following columns:
    I = 1:k-1; 
    J = k+1:n;
    scl = max(max(abs(E(:,k))), max(abs(A(:,k))));
    
    % Multiply the kth column of A with the basis in E:
    ex = innerprod(E(:,k), A(:,k));
    aex = abs(ex);
    
    % Adjust the sign of the kth column in E:
    if ( aex < eps*scl )
        s = 1; 
    else
        s = -sign(ex/aex);
    end
    E(:,k) = E(:,k) * s;
    
    % Compute the norm of the kth column of A:
    r = sqrt(innerprod(A(:,k), A(:,k)));
    R(k,k) = r;
    
    % Compute the reflection v:
    v = r*E(:,k) - A(:,k);
    % Make it more orthogonal:
    for j = I
        ev = innerprod(E(:,j), v);
        v = v - E(:,j)*ev;
    end
    
    % Normalize and store v:
    nv = sqrt(innerprod(v, v));
    if ( nv < tol*scl )
        v = E(:,k);
    else
        v = v / nv;
    end
    V(:,k) = v;
    
    % Subtract v from the remaining columns of A:
    for j = J
        av = innerprod(v, A(:,j));
        A(:,j) = A(:,j) - 2*v*av;
        rr = innerprod(E(:,k), A(:,j));
        A(:,j) = A(:,j) - E(:,k)*rr;
        R(k,j) = rr;
    end
    
end

% Form a discrete Q from the columns of V:
Q = E;
for k = n:-1:1
    for j = k:n
        vq = innerprod(V(:,k), Q(:,j));
        Q(:,j) = Q(:,j) - 2*V(:,k)*vq;
    end
end

% Compute the corresponding Chebysehv coefficients:
f.coeffs = f.chebpoly(Q);
% Trim the unneeded ones:
f.coeffs(1:n,:) = [];
% Comute new values:
f.values = f.chebpolyval(f.coeffs);

% Update the vscale:
f.vscale = max(abs(Q), [], 1);

% Additional output argument:
if ( nargout == 3 )
    if ( nargin == 2 && strcmp(flag, 'vector') )
        E = 1:m;
    else
        E = eye(m);
    end
end

end

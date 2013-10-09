function [Q, R] = abstractQR(A, E, myInnerProduct, myNorm, tol)
%ABSTRACTQR  Abstract implementation of householder method for QR factorisation.
%
% A is the input 'quasimatix' (or array-valued Chebfun, or matrix of values).
% E is a Legendre matrix with the same representation
% myInnerProduct(u, v) returns the L2 inner product of u and v.
% myNorm(v) returns the norm of v (or a suitable estimate).

% Pre-allocate the matrices R and V:
numCols = size(A, 2);
R = zeros(numCols);
V = A;                  % cols of V will store Householder vectors

if ( nargin < 5 )
    tol = eps;
end
if ( nargin < 4 )
    myNorm = @norm;
end

% Trefethen's stuff:
for k = 1:numCols

    % Indices of the previous and following columns:
    I = 1:k-1; 
    J = k+1:numCols;
    
    % Scale:
    scl = max(myNorm(E(:,k)), myNorm(A(:,k)));

    % Multiply the kth column of A with the basis in E:
    ex = myInnerProduct(E(:,k), A(:,k));
    aex = abs(ex);

    % Adjust the sign of the kth column in E:
    if ( aex < tol*scl )
        s = 1; 
    else
        s = -sign(ex/aex);
    end
    E(:,k) = E(:,k) * s;

    % Compute the norm of the kth column of A:
    r = sqrt(myInnerProduct(A(:,k), A(:,k)));
    R(k,k) = r;

    % Compute the reflection v:
    v = r*E(:,k) - A(:,k);
    % Make it more orthogonal:
    for i = I
        ev = myInnerProduct(E(:,i), v);
        v = v - E(:,i)*ev;
    end
    % Normalize:
    nv = sqrt(myInnerProduct(v, v));
    if ( nv < tol*scl )
       v = E(:,k); 
    else
       v = v / nv; 
    end
    % Store:
    V(:,k) = v;

    % Subtract v from the remaining columns of A:
    for j = J
        av = myInnerProduct(v, A(:,j));
        A(:,j) = A(:,j) - 2*v*av;
        rr = myInnerProduct(E(:,k), A(:,j));
        A(:,j) = A(:,j) - E(:,k)*rr;
        R(k,j) = rr;
    end

end

% Form Q from the columns of V:
Q = E;
for k = numCols:-1:1
    for j = k:numCols
        vq = myInnerProduct(V(:,k), Q(:,j));
        Q(:,j) = Q(:,j) - 2*(V(:,k)*vq);
    end
end

end


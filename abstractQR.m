function [Q, R] = abstractQR(A, E, myInnerProduct, myNorm, tol)
%ABSTRACTQR   Abstract implementation of Householder QR factorisation algorithm.
%   [Q, R] = ABSTRACTQR(A, E, MYINNERPRODUCT) computes a weighted QR
%   factorisation of A, where A is any "matrix-like" object that admits such a
%   decomposition, using an abstract implementation of the method in [1].  E is
%   a "matrix-like" object of the same type as A that functions as described in
%   [1], and MYINNERPRODUCT is a binary function which when given two "vectors"
%   of the type used to form the columns of A as arguments returns their L2
%   inner product.  MYINNERPRODUCT should be conjugate-linear in its first
%   argument.
%
%   [Q, R] = ABSTRACTQR(A, E, MYINNERPRODUCT, MYNORM) does the same but uses
%   MYNORM instead of NORM to estimate the sizes of the "vectors" for the
%   purposes of determining thresholds used in the algorithm.  This can save
%   computation time if NORM is expensive to compute but MYNORM provides a
%   cheaply computable estimate for the same result.  Note that neither NORM
%   nor MYNORM are used in the normalization of the "columns" of Q, which is
%   handled using the norm associated to MYINNERPRODUCT.
%
%   [Q, R] = ABSTRACTQR(A, E, MYINNERPRODUCT, MYNORM, TOL) does the same but
%   uses TOL instead of the default EPS when determining thresholds used in the
%   algorithm.
%
%   Example (QR of a random 3 x 3 matrix):
%     [Q, R] = abstractQR(randn(3), eye(3), @(u, v) u'*v);
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%  This function exists primarily to de-duplicate the Householder QR code,
%  which is used as the backbone of the QR implementations at all levels of the
%  system.  The main usage case is that A is a matrix of values on a grid used
%  as the foundation for some function representation technology.  In this
%  case, the matrix E will need to be a Legendre matrix for this representation
%  (e.g., a "Legendre-Chebyshev-Vandermonde" matrix if the grid consists of
%  Chebyshev points) so that the integrals underlying the inner product can be
%  computed exactly.  A can also be an array-valued CHEBFUN or CHEBTECH or a
%  quasimatrix, in which case E is typically an object of the same type whose
%  columns are Legendre polynomials.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        % Apply the Householder reflection:
        av = myInnerProduct(v, A(:,j));
        A(:,j) = A(:,j) - 2*v*av;

        % Compute other nonzero entries in the current row and store them:
        rr = myInnerProduct(E(:,k), A(:,j));
        R(k,j) = rr;

        % Subtract off projections onto the current vector E(:,k):
        A(:,j) = A(:,j) - E(:,k)*rr;
    end

end

% Form Q from the columns of V:
Q = E;
for k = numCols:-1:1
    for j = k:numCols
        % Apply the reflection again:
        vq = myInnerProduct(V(:,k), Q(:,j));
        Q(:,j) = Q(:,j) - 2*(V(:,k)*vq);
    end
end

end

function theta = subspace(A, B)
%SUBSPACE   Angle between subspaces.
%   SUBSPACE(A, B) finds the angle between two subspaces specified by the
%   columns (or rows, if A and B are transposed) of the CHEBFUN objects A and B.
%
%   If the angle is small, the two spaces are nearly linearly dependent.
%
%   References:
%   [1] A. Bjorck & G. Golub, Numerical methods for computing angles between 
%       linear subspaces, Math. Comp. 27 (1973), pp. 579-594.
%   [2] P.-A. Wedin, On angles between subspaces of a finite dimensional inner 
%       product space, in B. Kagstrom & A. Ruhe (Eds.), Matrix Pencils, Lecture 
%       Notes in Mathematics 973, Springer, 1983, pp. 263-285.
%   [3] A. V. Knyazev and M. E. Argentati, Principal Angles between Subspaces
%       in an A-Based Scalar Product: Algorithms and Perturbation Estimates.
%       SIAM Journal on Scientific Computing, 23 (2002), no. 6, 2009-2041.
%       http://epubs.siam.org:80/sam-bin/dbq/article/37733

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(A, 'chebfun') || ~isa(B, 'chebfun') )
    error('CHEBFUN:CHEBFUN:subspace:argin', ...
        'Both A and B must be column CHEBFUN.')
end
if ( ~domainCheck(A, B) )
    error('CHEBFUN:CHEBFUN:subspace:domain', ...
        'Domain mismatch.')
end
if ( A(1).isTransposed )
    if ( ~B(1).isTransposed )
        error('CHEBFUN:CHEBFUN:subspace:trans', ...
            'Dimension mismatch (transpose).');
    end
    A = A.';
    B = B.';
end

% Compute orthonormal bases of A and B:
A = orth(A);
B = orth(B);

% Compute the inner product:
C = innerProduct(A, B);

% Singular values of the inner product:
S = svd(C);

% Cosine of the angle: (see [3])
cosTheta = min(S); 

% Is the angle large? 
if ( cosTheta < 0.8 )
    theta = acos(cosTheta);
else 
    % If the angle is small, recompute using sine formulation:
    if ( size(A, 2) < size(B, 2) )
        sinTheta = norm(A - B*C');
    else
        sinTheta = norm(B - A*C);
    end
    theta = asin(sinTheta);
end

end

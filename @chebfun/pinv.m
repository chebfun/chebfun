function X = pinv(A, tol)
%PINV   Pseudoinverse of a column CHEBFUN.
%   X = PINV(A) produces a row CHEBFUN X so that A*X*A = A and X*A*X = X.
%
%   X = PINV(A, TOL) uses the tolerance TOL. The computation uses SVD(A) and any
%   singular value less than the tolerance TOL is treated as zero.
%
% See also SVD, RANK.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

if ( A(1).isTransposed ) 
    error('CHEBFUN:CHEBFUN:pinv:row', ...
        'PINV only defined for column CHEBFUN objects.')
end

% Compute the SVD:
[U, S, V] = svd(A, 0);
s = diag(S);

% Choose a tolerance if none is given:
if ( nargin == 1 )
	tol = max(length(A)*eps(max(s)), vscale(A)*eps);
end

% Compute the rank:
r = sum(s > tol);

if ( r == 0 )
	X = 0*A';
else
    U = extractColumns(U, 1:r);
    S = diag(ones(r,1)./s(1:r));
    V = V(:,1:r);
	X = V*S*U';
end

function Q = orth(A)
%ORTH   Array-valued CHEBFUN orthogonalization.
%   Q = ORTH(A) is an orthonormal basis for the range of the column CHEBFUN A.
%   That is, the columns of Q span the same space as the columns of A, Q'*Q = I,
%   and the number of columns of Q is the rank of A.
%
% See also NULL, SVD, RANK, QR.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( A(1).isTransposed ) 
	error('CHEBFUN:CHEBFUN:orth:row', ...
        'ORTH() only defined for column CHEBFUN objects.')
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

% Select these columns of U:
Q = extractColumns(U, 1:r);

end

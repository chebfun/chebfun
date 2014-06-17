function Z = null(A)
%NULL   Null space of an array-valued CHEBFUN.
%   Z = NULL(A) is an orthonormal basis for the null space of the column
%   CHEBFUN A.
%
% See also ORTH, SVD, RANK, QR.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( A(1).isTransposed )
	error('CHEBFUN:CHEBFUN:null:row', ...
        'NULL() only defined for column CHEBFUN objects.')
end

% Compute the SVD:
[U, S, V] = svd(A, 0);
s = diag(S);

% Choose a tolerance if none is given:
if ( nargin == 1 )
	tol = max(length(A)*eps(max(s)), vscale(A)*epslevel(A));
end

% Compute the rank:
r = sum(s > tol);

% Select null vectors from final columns of V:
dimV = size(V, 2);
Z = V(:,r+1:dimV);

end

function r = rank(A, tol)
%RANK   Rank of an array-valued CHEBFUN.
%   RANK(A) produces an estimate of the number of linearly independent columns
%   (or rows) of A.
%
%   RANK(A, TOL) is the number of singular values of A greater than TOL.
%
% See also SVD, QR.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Compute the singular values:
s = svd(A);

% Choose a tolerance if none is given:
if ( nargin == 1 )
	tol = max(length(A)*eps(max(s)), vscale(A)*eps);
end

% Compute the rank:
r = sum(s > tol);

end

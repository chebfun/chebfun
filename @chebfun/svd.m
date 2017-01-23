function [U, S, V] = svd(A, econ)
%SVD   Singular value decomposition of a CHEBFUN.
%   [U, S, V] = SVD(A) or SVD(A, 0), where A is an array-valued column CHEBFUN
%   with N columns, produces an N x N diagonal matrix S with nonnegative
%   diagonal elements in nonincreasing order, a column CHEBFUN U with N
%   orthonormal columns, and an N x N unitary matrix V such that A = U*S*V'.
%
%   If A is a row CHEBFUN with N rows, then U is a unitary matrix and V' is an
%   array-valued row CHEBFUN.
%
%   S = SVD(A) returns a vector containing the singular values of A.
%
%   The computation is carried out by orthogonalization operations following
%   Battles's 2005 thesis [1].
%
%   [1] Z. Battles, "Numerical Linear Algebra for Continuous Functions",
%   D. Phil. thesis, University of Oxford, 2005.
%
% See also QR, MRDIVIDE, RANK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( (nargin > 2) || ((nargin == 2) && (econ ~= 0)) )
    error('CHEBFUN:CHEBFUN:svd:twoArgs',...
          'Use svd(A) or svd(A, 0) for SVD of array-valued CHEBFUN.');
end

if ( A(1).isTransposed )    % A is a row CHEBFUN
    % Call CHEBFUN/QR():
    [Q, R] = qr(A');
    % Call discrete SVD():
    [V, S, U] = svd(R, 0);
    % Make V a CHEBFUN:
    V = Q*V;
    
else                        % A is a column CHEBFUN
    % Call CHEBFUN/QR():
    [Q, R] = qr(A);
    % Call discrete SVD():
    [U, S, V] = svd(R, 0);
    % Make U a CHEBFUN:
    U = Q*U;
    
end

% Output only the singular values:
if ( nargout < 2 )
    U = diag(S);
end

end

function varargout = lu(varargin)
%LU   LU factorization of a CHEBFUN2.
%
% [L, U] = LU( F ) returns two quasimatrices L and U of size inf by k and
% k by inf, respectively, where k is the rank of the CHEBFUN2 F.
% The quasimatrices L and U are "psychologically" lower and upper triangular.
% L is also unit lower triangular. This is computed by a continuous analogue of
% Gaussian elimination with complete pivoting.
%
% [L, U, P] = lu( F ) returns a k by 2 matrix P, containing the pivoting elements.
%
% [L, U, P, Q] = LU( F ) returns two vectors P and Q containing the x-values and
% y-values of the pivoting elements.
%
% [L, U, P] = LU(F, THRESH) returns the LU factorization of F that removes
% the tail of pivots below THRESH.
%
% For more information about the factorization:
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, submitted, 2014.
%
% See also CHOL, QR.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = lu@separableApprox(varargin{:});

end

function varargout = qr(varargin)
%QR Orthogonal-triangular decomposition of a CHEBFUN2.
%
% [Q, R] = QR( F ), where F is a separableApprox, produces an unitary column
% quasimatrix Q and a upper-triangular row quasimatrix R so that F = Q * R. This
% is computed by a continuous analogue of QR.
%
% [Q, R] = QR( F, 0 ) is the same as QR( F ).
%
% [Q, R, E] = QR( F ) and [Q, R, E] = QR( F, 'vector') produces a vector E that
% stores the pivoting locations.
%
% For more information about this decomposition:
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, Proc. Royal Soc. A., 2015.
%
% See also LU, and CHOL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = qr@separableApprox(varargin{:});

end

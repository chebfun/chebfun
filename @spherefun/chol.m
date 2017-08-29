function varargout = chol(varargin)
%CHOL    Cholesky factorization of a SPHEREFUN. 
%
% R = CHOL( F ), if F is a nonnegative definite SPHEREFUN then this returns
% an upper triangular quasimatrix so that R'*R is a decomposition of F. If
% F is not nonnegative definite then an error is thrown.
%
% L = CHOL(F, 'lower'), if F is a nonnegative definite SPHEREFUN then this
% produces a lower triangular quasimatrix so that L*L' is a decomposition
% of F. If F is not nonnegative definite then an error is thrown.
% 
% [R, p] = CHOL( F ), with two outputs never throwns an error message. If F
% is nonnegative definite then p is 0 and R is the same as above. If F is
% symmetric but negative definite or semidefinite then p is a positive
% integer such that R has p columns and R'*R is a rank p nonnegative
% definite SPHEREFUN that approximates F. This is particular useful when F
% is nonnegative definite, but rounding error have perturbed it to be
% semidefinite.
%
% [L, p] = CHOL(F, 'lower') same as above but the first argument is lower 
% triangular. 
%
% For more information about the factorization: 
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, Proc. Royal Soc. A., 2015. 
%
% See also SPHEREFUN/LU, and SPHEREFUN/QR. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = chol@separableApprox(varargin{:});

end

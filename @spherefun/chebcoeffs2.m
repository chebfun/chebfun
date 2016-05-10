function varargout = chebcoeffs2(varargin)
%CHEBCOEFFS2    Bivariate expansion coefficient
%   X = CHEBCOEFFS2(F) returns the matrix of bivariate coefficients such that
%       F= sum_{i=0}^{n-1} ( sum_{j=0}^{n-1} X(i+1,j+1) T_i(y) T_j(x) ).
%
%   [A, D, B] = CHEBCOEFFS2( f ) returns the same coefficients keeping them in
%   low rank form, i.e., X = A * D * B'.
%
%
% See also SPHEREFUN/PLOTCOEFFS2, SPHEREFUN/COEFFS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = chebcoeffs2@separableApprox(varargin{:});
end

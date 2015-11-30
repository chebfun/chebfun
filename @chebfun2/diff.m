function varargout = diff(varargin)
%DIFF   Derivative of a CHEBFUN2 object.
%   DIFF(F) is the derivative of F along the y direction.
%
%   DIFF(F, N) is the Nth derivative of F in the y direction.
%
%   DIFF(F, N, DIM) is the Nth derivative of F along the dimension DIM.
%     DIM = 1 (default) is the derivative in the y-direction.
%     DIM = 2 is the derivative in the x-direction.
%
%   DIFF(F, [NX NY]) is the partial derivative of NX of F in the first variable,
%   and NY of F in the second derivative. For example, DIFF(F,[1 2]) is
%   d^3F/dxd^2y.
%
% See also GRADIENT, SUM, PROD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = diff@separableApprox(varargin{:});

end

function varargout = del2(varargin)
%DEL2   Scaled Laplacian of a CHEBFUN2.
%   L = del2(F) is a numerical approximation of
%       del^2 F/4 = (d^2F/dx^2 + d^2F/dy^2)/4.
%
% See also LAPLACIAN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = del2@separableApprox(varargin{:});

end

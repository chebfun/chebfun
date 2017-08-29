function varargout = tan(varargin)
%TAN   Tangent of a SPHEREFUN.
%
% See also SPHEREFUN/SIN, SPHEREFUN/COS, SPHEREFUN/TANH

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = tan@separableApprox(varargin{:});

end

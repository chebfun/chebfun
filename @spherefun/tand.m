function varargout = tand(varargin)
%TAND  Tangent of a SPHEREFUN (in degrees)
%
% See also SPHEREFUN/TAN

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = tand@separableApprox(varargin{:});

end

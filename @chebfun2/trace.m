function varargout = trace(varargin)
% TRACE integral of a CHEBFUN2 along its diagonal
%
% TRACE(f) is the integral of function f(x,x).
%
% See also DIAG.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = trace@separableApprox(varargin{:});

end

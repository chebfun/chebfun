function varargout = trace(varargin)
% TRACE integral of a SPHEREFUN along its diagonal 
%
% TRACE(f) is the integral of function f(x,x) in spherical coordinates.
% 
% See also DIAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = trace@separableApprox(varargin{:});

end

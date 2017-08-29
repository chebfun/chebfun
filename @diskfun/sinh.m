function varargout = sinh(varargin)
%SINH   Hyperbolic sine of a DISKFUN.
%
% See also DISKFUN/SIN and DISKFUN/COSH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sinh@separableApprox(varargin{:});

end
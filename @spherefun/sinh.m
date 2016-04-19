function varargout = sinh(varargin)
%SINH   Hyperbolic sine of a SPHEREFUN.
%
% See also SPHEREFUN/SIN and SPHEREFUN/COSH

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = sinh@separableApprox(varargin{:});
end

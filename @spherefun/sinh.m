function varargout = sinh(varargin)
%SINH   Hyperbolic sine of a SPHEREFUN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = sinh@separableApprox(varargin{:});
end

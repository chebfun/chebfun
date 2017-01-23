function varargout = integral2(varargin)
%INTEGRAL2  Double integral of a SPHEREFUN over its domain.
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   SPHEREFUN.
%
% See also SPHEREFUN/INTEGRAL, SPHEREFUN/SUM2, SPHEREFUN/QUAD2D.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sum2(varargin{:});

end

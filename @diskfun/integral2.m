function varargout = integral2(varargin)
%INTEGRAL2  Double integral of a DISKFUN over its domain.
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   DISKFUN.
%
% See also DISKFUN/INTEGRAL, DISKFUN/SUM2, DISKFUN/QUAD2D.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sum2(varargin{:});

end

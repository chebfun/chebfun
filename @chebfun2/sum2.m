function varargout = sum2(varargin)
%SUM2   Double integral of a CHEBFUN2 over its domain.
%   I = SUM2(F) returns the double integral of a CHEBFUN2.
%
% See also INTEGRAL2, INTEGRAL, QUAD2D.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sum2@separableApprox(varargin{:});

end

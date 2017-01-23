function varargout = quad2d(varargin)
%QUAD2D  Compute definite integral of a CHEBFUN2.
%   I = QUAD2D( F ) returns the definite integral of a CHEBFUN2
%   over its domain of definition.
%
%   I = QUAD2D(F, a, b, c, d) returns the definite integral of a CHEBFUN2
%   over the rectangle [a,b] x [c,d].
%
% See also INTEGRAL2, SUM2, INTEGRAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = quad2d@separableApprox(varargin{:});

end

function varargout = quad2d(varargin)
%QUAD2D  Complete definite integral of CHEBFUN2.
%   I = QUAD2D( F ), returns the definite integral of a CHEBFUN2 integrated
%   over its domain of definition.
%
%   I = QUAD2D(F, a, b, c, d), returns the definite integral of a CHEBFUN2.
%   Integrated over the domangle [a b] x [c d].
%
% See also INTEGRAL2, SUM2, INTEGRAL.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = quad2d@separableApprox(varargin{:});

end

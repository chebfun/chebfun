function varargout = integral(varargin)
%INTEGRAL   Complete definite integral of CHEBFUN2.
%
%   I = INTEGRAL(F), returns the definite integral of a CHEBFUN2. Integrated
%   over its domain of definition.
%
%   I = INTEGRAL(F, g), returns the integral of a CHEBFUN2 along the curve
%   defined by the complex-valued CHEBFUN g.
%
% See also INTEGRAL2, SUM2, QUAD2D.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = integral@separableApprox(varargin{:});

end

function varargout = integral(varargin)
%INTEGRAL   Complete definite integral of SPHEREFUN2. 
%
%   I = INTEGRAL(F), returns the definite integral of a SPHEREFUN2. Integrated
%   over its domain of definition.
% 
%   I = INTEGRAL(F, g), returns the integral of a SPHEREFUN2 along the curve
%   defined by the complex-valued CHEBFUN g.
% 
% See also INTEGRAL2, SUM2, QUAD2D.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = integral@separableApprox(varargin{:});
end

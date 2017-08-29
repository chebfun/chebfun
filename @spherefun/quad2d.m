function varargout = quad2d(varargin)
%QUAD2D  Complete definite integral of SPHEREFUN. 
%   I = QUAD2D( F ), returns the definite integral of a SPHEREFUN integrated
%   over its domain of definition.
% 
%   I = QUAD2D(F, a, b, c, d), returns the definite integral of a SPHEREFUN.
%   Integrated over the domangle [a b] x [c d].
% 
% See also SPHEREFUN/INTEGRAL2, SPHEREFUN/SUM2, SPHEREFUN/INTEGRAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = quad2d@separableApprox(varargin{:});

end

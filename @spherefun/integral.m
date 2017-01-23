function varargout = integral(varargin)
%INTEGRAL   Complete definite integral of SPHEREFUN. 
%
%   I = INTEGRAL(F), returns the definite integral of a SPHEREFUN. Integrated
%   over its domain of definition.
% 
% See also SPHEREFUN/INTEGRAL2, SPHEREFUN/SUM2, SPHEREFUN/QUAD2D.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sum2(varargin{:});

end

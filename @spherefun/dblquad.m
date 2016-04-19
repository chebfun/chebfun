function varargout = dblquad(varargin)
%DBLQUAD   Complete definite integral of SPHEREFUN. 
%   I = DBLQUAD(F, a, b, c, d), returns the definite integral of a SPHEREFUN over
%   the region [a, b, c, d].
% 
%   This function is a wrapper for quad2d.
%
% See also SPHEREFUN/QUAD2D, SPHEREFUN/INTEGRAL2, SPHEREFUN/SUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = dblquad@separableApprox(varargin{:});

end

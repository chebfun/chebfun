function varargout = dblquad(varargin)
%DBLQUAD   Complete definite integral of CHEBFUN2.
%   I = DBLQUAD(F, a, b, c, d), returns the definite integral of a CHEBFUN2 over
%   the region [a, b, c, d].
%
%   This function is a wrapper for quad2d.
%
% See also QUAD2D, INTEGRAL2, SUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = dblquad@separableApprox(varargin{:});

end

function varargout = integral2(varargin)
%INTEGRAL2  Double integral of a CHEBFUN2 over its domain.
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   CHEBFUN2.
%
%   I = INTEGRAL2(F, [a b c d]) integrate F over the rectangle region [a b] x [c
%   d] provide this rectangle is in the domain of F.
%
%   I = INTEGRAL2(F, C) computes the volume under the surface F over the region
%   D with boundary C. C should be a complex-valued CHEBFUN that represents a
%   closed curve. This can be a very slow feature.
%
% See also INTEGRAL, SUM2, QUAD2D.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = integral2@separableApprox(varargin{:});

end

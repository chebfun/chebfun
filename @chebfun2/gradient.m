function varargout = gradient(varargin)
%GRADIENT  Gradient of a CHEBFUN2.
%   [FX FY] = GRADIENT(F) returns the gradient of the CHEBFUN2 F,
%   where FX is the derivative of F in the x direction and FY is the derivative
%   in the y direction. Both derivatives are returned as CHEBFUN2 objects.
%
%   G = GRADIENT(F) returns a CHEBFUN2V which represents
%
%            G = ( F_x ; F_y )
%
% See also GRAD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = gradient@separableApprox(varargin{:});

end

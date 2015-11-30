function varargout = grad(varargin)
%GRAD   Numerical gradient of a CHEBFUN2.
%   [FX FY] = GRAD(F) returns the numerical gradient of the CHEBFUN2 F, where FX
%   is the derivative of F in the x direction and FY is the derivative of F in
%   the y direction. Both derivatives are returned as CHEBFUN2 objects.
%
%   G = GRAD(F) returns a CHEBFUN2V which represents
%
%            G = (F_x ; F_y )
%
%  This command is shorthand for GRADIENT(F).
%
% See also GRADIENT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = grad@separableApprox(varargin{:});

end

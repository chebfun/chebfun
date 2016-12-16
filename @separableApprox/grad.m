function varargout = grad( f )
%GRAD   Numerical gradient of a SEPARABLEAPPROX.
%   [FX FY] = GRAD(F) returns the numerical gradient of the SEPARABLEAPPROX F, where FX
%   is the derivative of F in the x direction and FY is the derivative of F in
%   the y direction. Both derivatives are returned as SEPARABLEAPPROX objects.
%
%   G = GRAD(F) returns a CHEBFUN2V which represents
%
%            G = (F_x ; F_y )
%
%  This command is shorthand for GRADIENT(F).
%
% See also GRADIENT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call GRADIENT:
if ( nargout <= 1 )
    out = gradient( f );
    varargout = { out };
else
    [fx, fy] = gradient( f );
    varargout = {fx, fy};
end

end

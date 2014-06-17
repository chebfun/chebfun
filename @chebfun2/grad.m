function varargout = grad( f )
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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

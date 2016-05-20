function varargout = gradient(f)
%GRADIENT   Gradient of a CHEBFUN3.
%   [FX, FY, FZ] = GRAD(F) returns the gradient of the CHEBFUN3 object F, 
%   where
%   FX is the partial derivative of F in the first variable,
%   FY is the partial derivative of F in the second varaiable, and 
%   FZ is the partial derivative of F in the third variable.
%
%   G = GRAD(F) returns a CHEBFUN3V object which represents
%            G = (FX; FY; FZ).
%
% See also CHEBFUN3/GRAD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fx = diff(f, 1, 1);   % diff in the x-variable.
fy = diff(f, 1, 2);   % diff in the y-variable.
fz = diff(f, 1, 3);   % diff in the z-variable.

F = chebfun3v(fx, fy, fz);

if ( nargout <= 1 ) 
    varargout = {F};
else
    varargout = F.components; 
end

end
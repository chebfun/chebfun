function varargout = gradient(f)
%GRAD   Gradient of a CHEBFUN3.
%   [FX FY FZ] = GRAD(F) returns the gradient of the CHEBFUN3 
%   object F, where FX is the derivative of F in the first variable,
%   FY is the derivative of F in the second varaiable, and 
%   FZ is the derivative of F in the third variable.
%
%   G = GRAD(F) returns a CHEBFUN3V object which represents
%
%            G = (F_x ; F_y; F_z)
%
%   See also CHEBFUN3/GRAD.

fx = diff(f, 1, 1);   % diff in x-variable
fy = diff(f, 1, 2);   % diff in y-variable 
fz = diff(f, 1, 3);   % diff in z-variable 

F = chebfun3v(fx, fy, fz);

if ( nargout <= 1 ) 
    varargout = {F};
else
    varargout = F.components; 
end

end

function varargout = gradient(f)
%GRADIENT   Gradient of a CHEBFUN3.
%   [FX, FY, FZ] = GRAD(F) returns the gradient of the CHEBFUN3 object F, 
%   where
%   FX is the partial derivative of F in the first variable,
%   FY is the partial derivative of F in the second varaiable, and 
%   FZ is the partial derivative of F in the third variable.
%
%   G = GRAD(F) returns a CHEBFUN3V object which represents
%            G = [FX; FY; FZ].
%
% See also CHEBFUN3/GRAD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = chebfun3v(diffx(f), diffy(f), diffz(f));

if ( nargout <= 1 ) 
    varargout = {F};
else
    varargout = F.components; 
end

end
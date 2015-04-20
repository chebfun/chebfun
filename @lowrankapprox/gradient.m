function varargout = gradient( f ) 
%GRADIENT  Numerical gradient of a CHEBFUN2. 
%   [FX FY] = GRADIENT(F) returns the numerical gradient of the CHEBFUN2 F,
%   where FX is the derivative of F in the x direction and FY is the derivative
%   of F in the y direction. Both derivatives are returned as CHEBFUN2 objects.
%
%   G = GRADIENT(F) returns a CHEBFUN2V which represents
% 
%            G = ( F_x ; F_y )
%
% See also GRAD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fx = diff(f, 1, 2);   % diff in x-variable
fy = diff(f, 1, 1);   % diff in y-variable 

F = chebfun2v( fx, fy );

if ( nargout <= 1 ) 
    varargout = { F };
else
    varargout = F.components; 
end

end 

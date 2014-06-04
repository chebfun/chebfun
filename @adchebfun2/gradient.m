function F = gradient( f ) 
%GRADIENT  Numerical gradient of an ADCHEBFUN2. 
%   F = GRADIENT(F) returns the numerical gradient of the ADCHEBFUN2 F,
%   where FX is the derivative of F in the x direction and FY is the derivative
%   of F in the y direction.
%
% See also GRAD.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

fx = diff(f, 1, 2);   % diff in x-variable
fy = diff(f, 1, 1);   % diff in y-variable 

F = [fx, fy];

end 
function F = gradient( f ) 
%GRADIENT Numerical gradient of a chebfun2. 
% 
%  [FX FY]=GRADIENT(F) returns the numerical gradient of the chebfun2 F.
%  FX is the derivative of F in the x direction and
%  FY is the derivative of F in the y direction. Both derivatives
%  are returned as chebfun2 objects. 
%
%  G = GRADIENT(F) returns a chebfun2v which represents
% 
%            G = (F_x ; F_y )
%
% See also GRAD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

fx = diff(f, 1, 2);   % diff in x-variable
fy = diff(f, 1, 1);   % diff in y-variable 

F = chebfun2v( {fx, fy} );

end 
function [fx,fy,fz] = gradient( f ) 
%GRADIENT  Numerical surface gradient of a SPHEREFUN. 
%   [FX,FY,FZ] = GRADIENT(F) returns the numerical surface gradient of the
%   SPHEREFUN F, where FX, FY, and FZ are the x, y, and z componets of the 
%   the gradient.
%
% See also GRAD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fx = diff(f, 1);   % diff in x-variable
fy = diff(f, 2);   % diff in y-variable 
fz = diff(f, 3);   % diff in z-variable

% % Project onto the sphere which is done by (I - [x;y;z][x;y;z]')[fx;fy;fz],
% % where [x;y;z] is the unit normal vector (which is just the position
% % vector on the sphere).
% 
% % % Components of the normal vector
% x = spherefun(@(x,y,z) x, f.domain);
% y = spherefun(@(x,y,z) y, f.domain);
% z = spherefun(@(x,y,z) z, f.domain);
% 
% dotProd = x.*fx + y.*fy + z.*fz;
% xdot = x.*dotProd;
% ydot = y.*dotProd;
% zdot = z.*dotProd;
% 
% fx = fx - xdot;
% fy = fy - ydot;
% fz = fz - zdot;

end


function varargout = xyzsphere
%XYYSPHERE  Spherefun objects for x, y, and z on the surface of the sphere.
%   [X, Y, Z] = CHEB.XYZSPHERE returns spherefun objects for the functions 
%   @(x,y,z) x, @(x,y,z) y, and @(x,y,z) z defined on the surface of the
%   sphere. 
%
%   CHEB.XYZSPHERE is shorthand for the expressions 
%   SPHEREFUN(@(X,Y,Z) X), 
%   SPHEREFUN(@(X,Y,Z) Y), and 
%   SPHEREFUN(@(X,Y,Z) Z). 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = spherefun(@(x,y,z) x);
y = spherefun(@(x,y,z) y);
z = spherefun(@(x,y,z) z);

if ( nargout > 1 )
    
    % For the syntax [x, y, z] = cheb.xyzsphere:
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
    
else
    % Put 'x', 'y', and 'z' into the workspace:
    assignin('base', 'x', x)
    assignin('base', 'y', y)
    assignin('base', 'z', z)
    
    % And display the variables:
    display(x)
    display(y)
    display(z)
end
end
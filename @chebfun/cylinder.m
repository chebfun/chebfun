function varargout = cylinder( r )
%CYLINDER Generate cylinder. Surface revolution of a chebfun to form a chebfun2.
%
%   [X, Y, Z] = CYLINDER(R) forms the unit cylinder based revolving the 
%   function R about the z-axis. X, Y, and Z are chebfun2 objects such that
%   surf(X,Y,Z) displays the cylinder. 
%
%   F = CYLINDER(R) constructs the chebfun2v that represents the surface of
%   revolution. SURF(F) displays the cylinder.
%
%   Omitting output arguments causes the cylinder to be displayed with a SURF
%   plot.
%
% See also CHEBFUN2/SURF, CYLINDER, CHEBFUN2/SPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

d = [r.domain 0 2*pi];   % surface of revolution domain

f = chebfun2( @(x,y) feval(r, x), d );  
u = chebfun2( @(u,v) u, d );
v = chebfun2( @(u,v) v, d );

F = [ f.*sin(v) ; f.*cos(v) ; u ];    % surface of revolution.

if ( nargout == 0 )
    surf(F);  % plot
elseif ( nargout == 1 )
    % return chebfun2v object
    varargout = { F }; 
else
    % return as a parameterisation.
    varargout = { F(1), F(2), F(3) };   
end

end

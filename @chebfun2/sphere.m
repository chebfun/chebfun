function varargout = sphere( r )
%SPHERE   Generate a spherical surface. (Not necessarily a sphere!)
%   SPHERE(R), where R is a CHEBFUN2 on the domain [0, pi] x [0, 2*pi] plots the
%   "sphere" of radius R(th,phi).
%
%   [X, Y, Z]=SPHERE(R) returns X, Y, and Z as CHEBFUN2 objects such that
%   SURF(X,Y,Z) plots a sphere of radius R(th,phi).
% 
%   F = SPHERE(R) returns the CHEBFUN2V representing the sphere of radius R.
%   SURF(F) plots a sphere of radius R.
%
%   Omitting output arguments causes the SPHERE command to be displayed with a
%   SURF command and no outputs are returned.
%
% For the unit sphere: 
%   r = chebfun2(@(th, phi) 1+0*th, [0 pi 0 2*pi]);
%   F = sphere( r );   surf( F )
%
% For a sea shell:
%   r = chebfun2(@(th, phi) phi, [0 pi 0 2*pi]);
%   F = sphere( r ); surf( F )
% 
% See also CYLINDER, ELLIPSOID.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Sphere with radius r(th, phi).  
dom = [0 pi 0 2*pi]; 
th = chebfun2( @(th,phi) th, dom );
phi = chebfun2( @(th,phi) phi, dom );

x = r.*sin( th ).*cos( phi );
y = r.*sin( th ).*sin( phi );
z = r.*cos( th );

if ( nargout == 0 )
    surf( x, y, z), axis equal
elseif ( nargout == 1 )
    % Construct a CHEBFUN2V object:
    varargout = { [ x ; y ; z ] };  
else 
    varargout = { x, y, z };
end
    
end

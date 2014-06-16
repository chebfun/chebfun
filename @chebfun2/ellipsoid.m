function varargout = ellipsoid( a, b, c )
%ELLIPSOID  Generate an ellipsoid-like surface. (Not necessarily an ellipsoid!)
%   ELLIPSOID(A,B,C), where A, B, and C are CHEBFUN2 objects on the domain [0
%   pi]x[0 2*pi] plots the "ELLIPSOID" of semi axis lengths A(th,phi),
%   B(th,phi), and C(th,phi).
%
%   [X Y Z]=ELLIPSOID(A,B,C) returns X, Y, and Z as CHEBFUN2 objects such that
%   SURF(X,Y,Z) plots an ELLIPSOID of semi axis lengths A(th,phi), B(th,phi),
%   and C(th,phi).
% 
%   F = ELLIPSOID(A,B,C) returns the CHEBFUN2V representing the ELLIPSOID
%   SURF(F) plots the ELLIPSOID.
%
%   Omitting output arguments causes the ELLIPSOID command to be displayed with
%   a SURF command and no outputs are returned.
%
% For the ellipsoid: 
%   a = chebfun2(@(th,phi) 1+0*th,[0 pi 0 2*pi]);
%   F = ellipsoid(a,2*a,3*a); surf(F)
%
% For a badly shaped ship: 
%   a = chebfun2(@(th,phi) th,[0 pi 0 2*pi])
%   F = ellipsoid(a,2*a,cos(2*a)+.25); surf(F)
%
% See also SPHERE, CYLINDER.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    error('CHEBFUN:CHEBFUN2:ellipsoid:inputs', 'Incorrect input arguments.');
end

% ELLIPSOID with axis lengths a(th,phi), b(th,phi), c(th,phi).  
dom = [0 pi 0 2*pi]; 
th = chebfun2(@(th, phi) th, dom);
phi = chebfun2(@(th, phi) phi, dom);

x = a.*sin(th).*cos(phi);
y = b.*sin(th).*sin(phi);
z = c.*cos(th);

if ( nargout == 0 )
    surf(x, y, z)
    axis equal
elseif ( nargout == 1 )
    varargout = { [x ; y ; z] };
else
    varargout = { x, y, z };
end
    
end

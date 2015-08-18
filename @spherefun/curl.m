function [gx,gy,gz] = curl( f ) 
%CURL  Numerical surface curl of a scalar SPHEREFUN. 
%   [FX,FY,FZ] = CURL(F) returns the numerical surface curl of the scalar
%   SPHEREFUN F, where FX, FY, and FZ are the x, y, and z componets of the 
%   the surface curl.  The surface curl of a scalar function is defined to
%   be the curl of the f*n, where n is the normal vector to the sphere.
%   The result is a vector field that is tangent to the sphere.  This is
%   the generalization to the sphere of the definition of the curl of
%   scalar on a 2D plane.
%
%   TODO: Add support for computing the surface curl of a vector on the
%   sphere, i.e. [FX,FY,FZ] = CURL(GX,GY,GZ)
%
% See also DIVERGENCE, GRADIENT, VORTICITY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% On the sphere in Cartesian cooridnates the curl of a scalar is just
% cross(n,gradient(f)), where n is the normal vector (which is just
% [x,y,z]^T).

[fx,fy,fz] = gradient(f);

% Components of the normal vector
x = spherefun(@(x,y,z) x, f.domain);
y = spherefun(@(x,y,z) y, f.domain);
z = spherefun(@(x,y,z) z, f.domain);

% First component: -z*df/dy + y*df/dz
gx = -z.*fy + y.*fz;
% Second component: z*df/dx - x*df/dz
gy = z.*fx - x.*fz;
% Third component: -y*df/dx + x*df/dy
gz = -y.*fx + x.*fy;

end


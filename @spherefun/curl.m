function G = curl( f ) 
%CURL   Numerical surface curl of a scalar SPHEREFUN. 
%   G = CURL(F) returns a SPHEREFUNV G representing the numerical surface
%   curl of the scalar SPHEREFUN F.  The surface curl of a scalar
%   function is defined to be the curl of F*N, where n is the normal
%   vector to the sphere. The result is a vector field that is tangent to
%   the sphere.  This is the generalization to the sphere of the definition
%   of the curl of a scalar on a 2D plane.
%
% See also GRADIENT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Empty check.
if isempty( f )
    G = spherefunv;
    return;
end

% On the sphere in Cartesian cooridnates the curl of a scalar is just
% cross(n,gradient(f)), where n is the normal vector (which is just
% [x,y,z]^T).

F = gradient(f);
Fc = F.components;

% Components of the normal vector
x = spherefun(@(x,y,z) x, f.domain);
y = spherefun(@(x,y,z) y, f.domain);
z = spherefun(@(x,y,z) z, f.domain);

% First component: -z*df/dy + y*df/dz
gx = -z.*Fc{2} + y.*Fc{3};
% Second component: z*df/dx - x*df/dz
gy = z.*Fc{1} - x.*Fc{3};
% Third component: -y*df/dx + x*df/dy
gz = -y.*Fc{1} + x.*Fc{2};

G = spherefunv(gx,gy,gz);

end


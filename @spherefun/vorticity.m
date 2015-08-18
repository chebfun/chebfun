function V = vorticity( fx, fy, fz ) 
%VORTICITY  Numerical surface vorticity of a SPHEREFUN. 
%   V = VORTICITY(FX, FY, FZ) returns the numerical surface vorticity of the
%   vector field [FX FY FZ], where each component is a SPHEREFUN F and the
%   whole field is expressed in Cartesian coordinates. 
%   
%   Vorticity is the defined as the normal component of the surface curl of
%   the vector [FX FY FZ] to the sphere.  It is the generalization of the
%   standard 2D scalar vorticity for the surface of the sphere. FX, FY, and
%   FZ should be the components of the vector field expressed in cartesian
%   coordinates.
%
% See also VORT, DIVERGENCE, GRADIENT, CURL

%
% Compute the surface curl of the vector field.
%

% First component: d/dy(fz)-d/dz(fy)
gx = diff(fz, 2) - diff(fy, 3);
% Second component: d/dz(fx)-d/dx(fz)
gy = diff(fx, 3) - diff(fz, 1);
% Third component: d/dx(fy)-d/dy(fx)
gz = diff(fy, 1) - diff(fx, 2);

% Now dot with the normal vector on the sphere, which is just [x y z] 
% (i.e. the position vector on the sphere).
% Components of the normal vector
dom = fx.domain;
x = spherefun(@(x,y,z) x, dom);
y = spherefun(@(x,y,z) y, dom);
z = spherefun(@(x,y,z) z, dom);

V = x.*gx + y.*gy + z.*gz;

end
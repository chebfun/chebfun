% Wrapper for evaluating a function defined in Cartesian cooridnates as a
% function defined in spherical coordinates.
% f - smooth function in Cartesian coordinates on the sphere (x,y,z)
% lam,th - point to evaluate f in spherical coordinates.
% coord - Type of spherical coordinate system:
%         coord = 0 (co-latitude) --> -pi <= lam < pi, 0 <= th <= pi
%         coord = 1 (latitude)    --> -pi <= lam < pi, -pi/2 <= th <= pi/2
% Default is co-latitude.
function fdf = sphf2cartf(f,lam,th,coord)

% TODO: Create separate wrapper functions for latitude/co-latitude so that
% this if can be removed (thus improving performance.

if nargin == 3
    coord = 0;  % Default is to use co-latitude.
end

if coord == 0
    % Latitude: 0 <= theta < pi
    x = cos(lam).*sin(th);
    y = sin(lam).*sin(th);
    z = cos(th);
elseif coord == 1
    % Latitude: -pi/2 <= theta < pi/2
    x = cos(lam).*cos(th);
    y = sin(lam).*cos(th);
    z = sin(th);
else
    error('SPHEREFUN:sphf2cartf:CoordSysUnknown: Unknown coordinate system for the sphere.');
end

fdf = f(x,y,z);

end
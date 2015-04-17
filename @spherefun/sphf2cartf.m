% Wrapper for evaluating a function defined in Cartesian cooridnates as a
% function defined in spherical coordinates.
% f - smooth function in Cartesian coordinates on the sphere (x,y,z)
% lam,th - point to evaluate f in spherical coordinates.
% coord - Type of spherical coordinate system:
%         coord = 0 (default) --->  -pi <= lam < pi, -pi/2 <= th <= pi/2
%         coord = 1 (co-latitude) --> 0 <= lam < 2*pi, 0 <= th <= pi
function fdf = sphf2cartf(f,lam,th,coord)

if nargin == 3
    coord = 0;  % Default is to not use co-latitude.
end

if coord == 1
    % Latitude: 0 <= theta < pi
    x = cos(lam).*sin(th);
    y = sin(lam).*sin(th);
    z = cos(th);
else
    % Latitude: -pi/2 <= theta < pi/2
    x = cos(lam).*cos(th);
    y = sin(lam).*cos(th);
    z = sin(th);
end


fdf = f(x,y,z);

end
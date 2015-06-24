% Wrapper for evaluating a function defined in Cartesian cooridnates as a
% function defined in spherical coordinates.
% f - smooth function in Cartesian coordinates on the sphere (x,y,z)
% lam,th - point to evaluate f in spherical coordinates.
% coord - Type of spherical coordinate system:
%         coord = 0 (co-latitude) --> -pi <= lam < pi, 0 <= th <= pi
%         coord = 1 (latitude)    --> -pi <= lam < pi, -pi/2 <= th <= pi/2
% Default is co-latitude.
function fdf = pol2cartf(f,th,r)

    x = r.*cos(th);
    y = r.*sin(th);
    

fdf = f(x,y);

end
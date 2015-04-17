% Wrapper for creating a "double Fourier" version of a function f defined
% on the sphere.
% f - smooth function in spherical coordinates on the sphere (lon,lat)
% lam - longitudinal (or azimuthal) point to evaluate f.
% th - latitudinal (or elevation) point to evaluate f.
% coord - Type of spherical coordinate system:
%         coord = 0 (default) --->  -pi <= lam < pi, -pi/2 <= th <= pi/2
%         coord = 1 (co-latitude) --> 0 <= lam < 2*pi, 0 <= th <= pi
function fdf = sph2torus(f,lam,th,coord)

if nargin == 3
    coord = 0;  % Default is to not use co-latitude.
end

fdf = real(f(lam,th));

if coord == 1
    id = th - pi > 100*eps;
    if ~isempty(id) && any(id(:))
        fdf(id) = f(mod(lam(id)-pi,2*pi),2*pi-th(id));
    end
else
    id = th-pi/2 > 100*eps;
    if ~isempty(id) && any(id(:))
        fdf(id) = f(mod(lam(id),2*pi)-pi,pi-th(id));
    end
end

end
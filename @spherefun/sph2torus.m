function fdf = sph2torus(f, lam, th, coord)
%SPH2TORUS  Wrapper for creating a "double Fourier" version of a function f
%   defined on the sphere.
%
%   f - smooth function in spherical coordinates on the sphere (lam, th)
%   lam - longitudinal (or azimuthal) point to evaluate f.
%   th  - latitudinal (or elevation) point to evaluate f.
%   coord - Type of spherical coordinate system:
%             coord = 0 (co-latitude) --> -pi <= lam < pi, 0 <= th <= pi
%             coord = 1 (latitude)    --> -pi <= lam < pi, -pi/2 <= th <= pi/2
%   Default is co-latitude

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% TODO: Create separate wrapper functions for latitude/co-latitude so that
% this if can be removed (thus improving performance.

if ( nargin == 3 )
    coord = 0;  % Default is to use co-latitude.
end

fdf = real(f(lam, th));

if ( coord == 0 )
    id = th - pi > 100*eps;
    if ( ~isempty(id) && any(id(:)) )
        fdf(id) = f(mod(lam(id)-pi,2*pi),2*pi-th(id));
    end
elseif ( coord == 1 )
    id = th - pi/2 > 100*eps;
    if ( (~isempty(id)) && (any(id(:))) )
        fdf(id) = f(mod(lam(id), 2*pi)-pi, pi-th(id));
    end
else
    error('SPHEREFUN:sph2cart:coordSysUnknown',...
        ['Unknown coordinate system for the sphere.']);
end

end
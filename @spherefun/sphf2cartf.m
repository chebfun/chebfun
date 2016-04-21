function g = sphf2cartf(f, lam, th, coord)
%SPHF2CARTF    Wrapper for evaluating a function defined in Cartesian 
% 
%  G = SPHF2CARTF(F, LAM, TH) evaluates the function handle F = F(x,y,z)
%  at x = cos(lam).*sin(th), y = sin(lam).*sin(th), z = cos(th). This is
%  the co-latitude spherical coordinate system.
%
%  G = SPHF2CARTF(F, LAM, TH, 0) same as SPHF2CARTF(F, LAM, TH).
%
%  G = SPHF2CARTF(F, LAM, TH, 1) evaluates F = F(x,y,z) at 
%  at x = cos(lam).*cos(th), y = sin(lam).*cos(th), z = sin(th). This is
%  the latitude spherical coordinate system.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% For Developers, recall that: 
%   coord - Type of spherical coordinate system:
%          coord = 0 (co-latitude) --> -pi <= lam < pi, 0 <= th <= pi
%          coord = 1 (latitude)    --> -pi <= lam < pi, -pi/2 <= th <= pi/2
% (Default above is co-latitude.)

% TODO: Create separate wrapper functions for latitude/co-latitude so that
% this if can be removed (thus improving performance.

if ( nargin == 3 )
    coord = 0;  % Default is to use co-latitude.
end

if ( coord == 0 )
    % Latitude: 0 <= theta < pi
    x = cos(lam).*sin(th);
    y = sin(lam).*sin(th);
    z = cos(th);
elseif ( coord == 1 )
    % Latitude: -pi/2 <= theta < pi/2
    x = cos(lam).*cos(th);
    y = sin(lam).*cos(th);
    z = sin(th);
else
    error('SPHEREFUN:sphf2cartf:CoordSysUnknown', ['Unknown coordinate '...
        'system for the sphere.']);
end

g = feval(f, x, y, z);

end
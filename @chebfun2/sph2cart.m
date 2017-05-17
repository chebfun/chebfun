function [x, y, z] = sph2cart( th, phi, r )
%SPH2CART   Transform spherical to Cartesian coordinates for CHEBFUN2 objects.
%   [X, Y, Z] = SPH2CART(TH, PHI, R) transforms corresponding elements of data
%   stored in spherical coordinates (azimuth TH, elevation PHI, radius R) to
%   Cartesian coordinates X,Y,Z.  The arrays TH, PHI, and R must be the same
%   size (or any of them can be scalar).  TH and PHI must be in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
% See also POL2CART.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers. 
% See http://www.chebfun.org/ for Chebfun2 information. 

% Spherical coordinates: 
z = r .* sin( phi );
rcos = r .* cos( phi );
x = rcos .* cos( th );
y = rcos .* sin( th );

end

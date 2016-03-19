function Y = sphharm(l, m, coord)
%SPHHARM    Normalized, real-valued, spherical harmonic of degree L, order 
%   M at a given set of locations on the sphere.
%
%   Y = sphHarm(L, M, COORD) returns the degree L, order M normalized
%   spherical harmonic on the sphere expressed in longitude-latitude
%   coordinates (or azimuthal-elevation).  Here
%        -pi <= lam <= pi   is the longitude (azimuthal) coordinate.
%   The flag `coord` determines whether co-latitude or latitude is to be
%   used:
%     coord = 0 ---> co-latitude 0 <= th <= pi (default)
%     coord = 1 ---> latitude -pi/2 <= th <= pi/2 

% The degree l must be greater than or equal to the magnitude of the order
% m
if ( l < abs(m) )
    error('SPHEREFUN:sphHarm', ['The degree of the spherical harmonic '...
        'must be greater than or equal to the magnitude of the order']);
end

if ( abs(l) > 75 )
    error('SPHEREFUN:sphHarm','Only spherical harmonics of degree <= 75.');
end

if ( nargin <= 2 )
    coord = 0;
end

if ( coord == 1 )
    dom = [-pi pi -pi/2 pi/2];
else
    dom = [-pi pi 0 pi];
end

Y = spherefun(@(lam, th) mySphHarm(l, m, lam, th, coord), dom);

end

function Y = mySphHarm(l, m, lam, th, coord)
%SPHHARM    Normalized, real-valued, spherical harmonic of degree
%   L, order M at a given set of locations on the sphere.
%
%   Y = sphHarm(L, M, LAM, TH) returns the degree L, order M normalized
%   spherical harmonic at the points (LAM, TH) on the sphere expressed in
%   longitude-latitude coordinates or azimuthal-elevation (lam,th).  Here
%        -pi <= lam <= pi   is the longitude (azimuthal) coordinate
%
%   The flag `coord` determines wether co-latitude or latitude is to be
%   used:
%     coord = 0 ---> co-latitude 0 <= th <= pi (default)
%     coord = 1 ---> latitude -pi/2 <= th <= pi/2 

%   TODO: this implementation is not stable or fast for large (l, m).  Alex,
%   can you come up with something better?

if ( nargin <= 4 )
    coord = 0;
end

if ( coord == 1 )
    z = sin(th);  % Latitude
else
    z = cos(th);  % Co-latitude
end
    
% Flatten and transpose th and lam so they work with the legendre function
sz = size(z); 
z = z(:)'; 
lam = lam(:)';

% Normalization
% a = sqrt((2*l+1)/2/pi*factorial(l-abs(m))/factorial(l+abs(m))*(2-double(m==0)));
kk = l-abs(m)+1:l+abs(m);
aa = exp(-sum(log(kk)));
a = sqrt((2*l+1)/2/pi * aa * (2 - double(m==0)));
Y = legendre(l, z);

% Get the right associated legendre function
Y = squeeze(Y(abs(m)+1, :, :));

% Determine if the cos or sin term should be added.
pos = abs(max(0, sign(m+1)));

% Compute the spherical harmonic
Y = (pos*cos(m*lam) + (1-pos)*sin(m*lam)).*(a*Y);

% Reshape so it is the same size as the th and lam that were passed in.
Y = reshape(Y,sz);

end
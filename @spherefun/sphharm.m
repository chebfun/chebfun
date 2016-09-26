function Y = sphharm(l, m)
%SPHHARM   Normalized, real-valued, spherical harmonic of degree L, order 
%   M at a given set of locations on the sphere.
%
%   Y = sphHarm(L, M) returns the degree L, order M normalized
%   spherical harmonic on the sphere expressed in longitude-latitude
%   coordinates (or azimuthal-elevation).  Here
%        -pi <= lam <= pi   is the longitude (azimuthal) coordinate, and
%          0 <= th  <= pi   is the latitude (elevation) coordinate.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Add support for constructing spherical harmonics with the latitude
% coordinate being -pi/2 <= th <= pi/2.

% The degree l must be greater than or equal to the magnitude of the order
% m:
if ( l < abs(m) )
    error('SPHEREFUN:sphHarm', ['The degree of the spherical harmonic '...
        'must be greater than or equal to the magnitude of the order']);
end

if ( abs(l) > 75 )
    error('SPHEREFUN:sphHarm','Only spherical harmonics of degree <= 75.');
end

% Developer note: once support for different longitude and latitude domains
% is added to spherefun, the following code below with coord = 0 or 
% coord = 1 can be used.  This coord parameter should be an optional 
% argument to the sphharm method.  For now we just always set coord to 0
% since this cooresponds to the only domain currently supported in
% spherefun.
coord = 0;

if ( coord == 1 )
    dom = [-pi pi -pi/2 pi/2];
else
    dom = [-pi pi 0 pi];
end

% We could construct a spherical harmonic function by calling the 
% Spherefun constructor using 
% Y = spherefun(@(lam, th) mySphHarm(l, m, lam, th, coord), dom);
% This would adaptively sample the function and determine "optimal" Fourier
% degrees to represent the spherical harmonic.  However, based on the 
% properties of the spherical harmonics, we already know ahead of time the
% correct number of samples to use.  Each function in longtidue (lambda) is
% a Fourier mode of degree m.  So, we can resolve this with 2|m|+2 samples.
% Each function in latitude (theta) is a associated Legendre function of
% theta. Specifically, it is a polynomial of cos(theta) or sin(theta)
% if coord=0 or coord=1, respectively. The highest degree this
% polynomial can be is l when when m=0.  So, we can resolve this using 2l+2
% samples.  

% Construct a matrix of values at the latitude-longitude grid
[ll,tt] = meshgrid(trigpts( max( 2*abs(m)+2, 4 ), dom(1:2) ),...
                   linspace( dom(3), dom(4), 2*l+2 ) );
Y = spherefun( mySphHarm(l,m,ll,tt,coord), dom );
% Simplify to get the most compressed representation.
Y = simplify( Y );

end

function Y = mySphHarm(l, m, lam, th, coord)
%MYSPHHARM   Main subroutine for computing the spherical harmonics.
%
% [TODO]: This implementation is not stable or fast for large (l, m) for the
% following reasons:
%
% 1. Stability: it uses the normalized spherical harmonics and the
% normalization factors are unstable to compute for large l and m.  We
% should switch to unnormalized spherical harmonics.
% 
% 2. Efficiency: the code just uses matlab's `legendre` function for 
% computing the associated legendre polynomials and this function is dead
% slow. Instead a recursion formula should be used.

if ( nargin <= 4 )
    coord = 0;
end

if ( coord == 1 )
    z = sin(th);  % Latitude
else
    z = cos(th);  % Co-latitude
end
    
% Flatten and transpose th and lam so they work with the legendre function:
sz = size(z); 
z = z(:)'; 
lam = lam(:)';

% Normalization:
kk = l-abs(m)+1:l+abs(m);
aa = exp(-sum(log(kk)));
a = sqrt((2*l + 1)/4/pi * aa * (2 - double(m==0)));
Y = legendre(l, z);

% Get the right associated legendre function:
Y = squeeze(Y(abs(m)+1, :, :));

% Determine if the cos or sin term should be added:
pos = abs(max(0, sign(m+1)));

% Compute the spherical harmonic:
Y = (pos*cos(m*lam) + (1-pos)*sin(m*lam)).*(a*Y);

% Reshape so it is the same size as the th and lam that were passed in.
Y = reshape(Y, sz);

end
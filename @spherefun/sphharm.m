function Y = sphharm(l, m)
%SPHHARM   Real-valued, spherical harmonic of degree L, order M.
%
%   Y = SPHHARM(L, M) returns the degree L, order M real-valued 
%   spherical harmonic on the sphere.  Y is normalized so that its two-norm
%   over the sphere is 1.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Add support for constructing spherical harmonics with the latitude
% coordinate being -pi/2 <= th <= pi/2.

% The degree l must be greater than or equal to the magnitude of the order
% m:
if ( l < abs(m) )
    error('CHEBFUN:SPHEREFUN:sphHarm', 'Degree must be >= order for spherical harmonic ');
end

if ( abs(l) > 120 )
    warning('CHEBFUN:SPHEREFUN:sphHarm:maxDegree','Results may be inaccurate for degrees larger than approximately 120.');
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
ll = trigpts( max( 2*abs(m)+2, 4 ), dom(1:2) );
tt = linspace( dom(3), dom(4), 2*l+2 );
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
% normalization factors are unstable to compute for large l and m.
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

% Make lam and z row vectors for what follows.
z = z.';
% Make lam row
lam = lam(:).';

% Get the normalized associated Legendre function.
Y = (-1)^m/sqrt((1+double(m==0))*pi)*legendre(l, z, 'norm');

% Get the right associated legendre function:
Y = squeeze(Y(abs(m)+1, :, :)).';

% Determine if the cos or sin term should be added:
pos = abs(max(0, sign(m+1)));

% Compute the spherical harmonic:
Y = Y*(pos*cos(m*lam) + (1-pos)*sin(abs(m)*lam));

end

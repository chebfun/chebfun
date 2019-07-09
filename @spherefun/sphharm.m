function Y = sphharm(l, m)
%SPHHARM   Real-valued, spherical harmonic of degree L, order M.
%
%   Y = SPHHARM(L, M) returns the degree L, order M real-valued 
%   spherical harmonic on the sphere.  Y is normalized so that its two-norm
%   over the sphere is 1. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( l < abs(m) )
    error('CHEBFUN:SPHEREFUN:sphHarm', 'Degree must be >= order for spherical harmonic ');
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
ll = trigpts( 2*abs(m)+2, dom(1:2) );
tt = linspace( dom(3), dom(4), 2*l+2 );
Y = spherefun( mySphHarm(l, m, ll, tt), dom );

end


function Y = mySphHarm(l_max, m_max, lam, th)
%MYSPHHARM   Main subroutine for computing the spherical harmonics.
%
% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%
% Below is the implementation the Modified Forward Column (MFC) method 
% described in the Holmes and Featherstones paper (2002)
% It computes P^m_n/u^m by a stable recurrence to avoid numerical errors
% near the poles

abs_m_max = abs( m_max );

% Make lam row
lam = lam.';
% Make th vector
th = th.';

% Precompute vector cos(th)
p = length(th);
CosTh = cos(th);

%% Compute P_m_max^m_max / u^m_max

% Initialize P^0_0 / u^0
Pold = ones(p, 1);

% Compute P^m_m/u^m
for m = 1:abs_m_max
    % Compute Pm^m with Pm-1^m-1
    Pold = sqrt( (2*m+1)/(2*m-(m==1)) ) * Pold;
end

% Initialize the recurrence (Pm^m-1 does not exist)
Poldold = zeros(p, 1);

%% Compute P^m_l / u^m with the recurrence formula, m_max+1 <= l <= l_max
for l = abs_m_max+1:l_max
    anm = sqrt( (4*l^2-1)/((l-abs_m_max)*(l+abs_m_max)) );
    bnm = sqrt( (2*l+1)*(l+abs_m_max-1)*(l-abs_m_max-1)/((l-abs_m_max)*(l+abs_m_max)*(2*l-3)) );
    % Compute the normalized associated legendre polynomial P^m_l/u^m
    Pl = anm*CosTh.*Pold - bnm*Poldold;

    % Update the polynomials for the recurrence
    Poldold = Pold;
    Pold = Pl;
end

% Normalize the polynomial and recover associated Legendre polynomials
Pold = (-1)^abs_m_max*sin(th).^abs_m_max.*Pold/sqrt(4*pi);

% Determine if the cos or sin term should be added:
pos = abs( max(0, sign(m_max+1)) );

% Compute the spherical harmonic:
Y = Pold*( pos*cos(m_max*lam) + (1-pos)*sin(abs(m_max)*lam) );

end

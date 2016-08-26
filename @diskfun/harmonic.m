function Y = harmonic(L, m, type)
%HARMONIC   Normalized, real-valued, fixed-height cylindrical harmonic.
%
%   Y = HARMONIC(L, M) returns the cylindrical harmonic V_L^M(t,r) with 
%   homogeneous Dirichlet boundary conditions. Here,
%        -pi <= t <= pi  is the angular coordinate, and
%          0 <= r  <= 1  is the radial coordinate.
%
%   Y = HARMONIC(L, M, 'dirichlet') is the same as Y = HARMONIC(L, M).
%
%   Y = HARMONIC(L, M,'neumann') returns the cylindrical harmonic V_L^M(t,r)
%   with Neumann boundary conditions.
%
%   The cylindrical harmonics are an orthogonal basis of functions in
%   cylindrical coordinates. If the height variable is set as a fixed
%   value, then cylindrical harmonics are the eigenfunctions of the
%   Laplacian on the disk. For Dirichlet boundary conditions, they are of
%   the form V_L^M(t, r) = A*exp(iLt).*J_L(a_Mr), where J_L is the Lth
%   j-Bessel function, a_M is the Mth positive zero of the function, and A
%   is a normalization factor.
%
% See also SPHEREFUN/SPHHARM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 || isempty(type) )
    type = 'dirichlet';
end

% Create a diskfun object for cylindrical harmonic:
Y = diskfun(@(theta, r) cylindricalharmonic(theta, r, L, m, type), 'polar');

end

function Z = cylindricalharmonic(theta, r, L, m, type)
%CYLINDRICALHARMONIC   Construct function handle for cylindrical harmonic.
%
%   Z = CYLINDRICALHARMONIC( THETA, R, L, M, TYPE ), construct a function
%   handle in coordinates (THETA, R) for the (L,M) cylindrical harmonic of
%   type TYPE. If TYPE = 'dirichlet', then the harmonic has homogeneous 
%   Dirichlet data; otherwise, TYPE = 'neumann', then the harmonic has 
%   homogeneous neumann data at the boundary. 

% Calculate the mth positive zero of the Lth Bessel function
Lsign = sign( L+1 );
L = abs( L );

if ( strcmpi(type, 'dirichlet') )
    % We would like to do the following (which is much faster). However,
    % that algorithm is not particularly accurate. 
    %  Jzero = besselroots(L, m);
    %  Jzero = Jzero( m );
    
    % Instead, we currently use chebfun/roots:
    Jzero = roots(chebfun(@(x) besselj(L,x), [sqrt((3/4)^2*pi^2+L^2) (m+L/2)*pi]));
    Jzero = Jzero(m);
    z = besselj(L, r*Jzero);

    % Choose sin or cos based on L value
    pos = abs( max(0, Lsign) );
    Z = ((pos)*cos(L*theta) + (1-pos)*sin(L*theta)).*z;
    
    % Normalize
    k = double( L==0 );
    Z = sqrt(2)/ (sqrt((1+k)*pi)*abs(besselj(L+1, Jzero)))*Z;
elseif ( strcmpi(type, 'neumann') )    
    if ( L==0 && m==1 )
        % Case where eigenvalue is zero; constant mode
        Z = 1 + r*0;  
    else
        % The x-term is not present when L=0
        Jzero = roots(chebfun(@(x) L*besselj(L,x)-x.^(min(1,L))...
            .*besselj(L+1, x), [L (m+L/2)*pi])); 
        Jzero = Jzero(m); % When L=0 this is still the correct index
        z = besselj(L, r*Jzero);
        % Choose sin or cos based on L value
        pos = abs(max(0,Lsign));
        Z = ((pos)*cos(L*theta) +(1-pos)*sin(L*theta)).*z;
        % Normalize
        k = double(L==0);
        Z = sqrt(2)/ (sqrt( (1-L^2/Jzero^2)*(1+k)*pi)*abs(besselj(L, Jzero)))*Z;
    end   
else    
    error('DISKFUN:HARMONIC:TYPE', 'Unrecognized Bessel type.')
end

end
function f = rotate(f, phi, theta, psi, method)
%ROTATE   Rotates a SPHEREFUN using Euler angles
%   Y = ROTATE(F, PHI, THETA, PSI) rotates F using Euler angles phi, theta, 
%   and psi with the ZXZ convention:  Rotate first about the z-axis by an
%   angle phi, then about the (orginal) x-axis by an angle 0<=theta<=pi, 
%   then about new z-axis by an angle psi. 
%
%   Y = ROTATE(F, PHI, THETA, PSI, METHOD) is the same as 
%   Y = ROTATE(F, PHI, THETA, PSI) except it uses the algorithm METHOD. The
%   options for METHOD are: 
%               - 'feval': Uses Horner's scheme for evaluation,
%               - 'nufft': Uses the 2D nonuniform FFT for evaluation.
%
% See also spherefun.fastSphereEval, diskfun.rotate.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    return
elseif ( nargin == 2 )
    theta = 0;
    psi = 0;
elseif ( nargin == 3 )
    psi = 0;
end

[m, n] = length( f );
if ( nargin == 4 )
    % Use the fast transform, unless the user tells us otherwise.
    method = 'nufft'; 
end

% f is bandlimited of degree (n,m) so its rotation will be bandlimited will
% limit at most max(n,m). This follows from spherical harmonic theory as
% the rotation of a spherical harmonic of degree l is of degree l, only the
% order will change. The degree l is given by the bandlimit m, while the
% order is given by the bandlimit of n.
n = max(m,n);             % Using this sampling we should exactly recover the rotated f.
m = n + mod(n,2);         % Number of columns must be even.
% Set the number of rows in the sampled grid equal to n/2+2 so that the 
% doubled up grid will have n rows (n/2+2 because the pole is included in
% the sampled grid, and the doubled up grid does not contain -pi).
n = ceil(n/2)+2;

% Sampling grid.
[lam,th] = meshgrid(trigpts(m,[-pi pi]),linspace(0,pi,n));
lam(1,:) = 0;
lam(n,:) = 0;
x = cos(lam).*sin(th);
y = sin(lam).*sin(th);
z = cos(th);

% Rotation built up from zxz Euler angle rotations
D = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
C = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
B = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
R = B*C*D;

% Rotate the sampling grid
u = R(1,1)*x + R(1,2)*y + R(1,3)*z;
v = R(2,1)*x + R(2,2)*y + R(2,3)*z;
w = R(3,1)*x + R(3,2)*y + R(3,3)*z;

% Get the spherical coordinates of the rotated grid
[lam, th] = cart2sph(u, v, w);
th = pi/2-th;  % Adjust elevation angle since matlab uses latitude.
 
if ( strcmpi(method, 'nufft') )     % NUFFT evaluation using 2D NUFFT
    
    g = real( fastSphereEval(f, lam, th) );
    f = spherefun( g );
    
elseif ( strcmpi(method, 'feval') ) % FEVAL evalutation using Horner
   
    g = feval(f, lam, th);
    f = spherefun( g );
    
else
    
    error( 'Unrecognized algorithm.' )

end

% Simplify the result
f = simplify( f );

end
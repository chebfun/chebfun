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
% See also spherefun.FastSphereEval, diskfun.rotate.

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

[m, n] = length(f);
if ( nargin == 4 )
    % Use the fast transform only when the number lengths of both columns
    % are greater than 64 and the rank of f is greater than 2.
    if ( min(m,n) > 64 && numel(cdr(f)) > 2) 
        method = 'nufft';
    else
        method = 'feval';
    end
end

% % Slowest of all rotates (never used):
% if type == 3
%     fr = @(lam,th) rotateColat(f, lam, th, alpha, theta);
%     f = spherefun(fr, f.domain);
%     return
% end

% f is bandlimited of degree (m,n) so its rotation will be bandlimited will
% limit at most max(m,n). This follows from spherical harmonic theory as th
% rotation of a spherical harmonic of degree l is of degree l, only the order
% will change. The degree l is given by the bandlimit n, while the order is
% given by the bandlimit of m. Numerically, we need 2*max(m,n). Why?
n = 2*max(m,n);             % Using this sampling we should exactly recover the rotated f.
n = n + mod(n,2);           % Number of columns must be even.
m = n/2 + 1 - mod(n/2,2);   % Number of rows must be odd

% Sampling grid.
[lam,th] = meshgrid(trigpts(n,[-pi pi]),linspace(0,pi,m));
lam(1,:) = 0;
lam(m,:) = 0;
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
    
elseif ( strcmpi(method, 'nufft2') ) % FEVAL evalutation using Horner
   
    g = real( fastSphereEval2(f, lam, th) );
    f = spherefun( g );

elseif ( strcmpi(method, 'feval') ) % FEVAL evalutation using Horner
   
    g = feval(f, lam, th);
    f = spherefun( g );
    
else
    
    error('Unrecognized algorithm.')

end

end

% % Slowest of all rotates (never used):
% function y = rotateColat(f, lam, th, phi, theta, psi)
%     x = cos(lam).*sin(th);
%     y = sin(lam).*sin(th);
%     z = cos(th);
%     D = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
%     C = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
%     B = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
%     R = B*C*D;
% 
%     u = R(1,1)*x + R(1,2)*y + R(1,3)*z;
%     v = R(2,1)*x + R(2,2)*y + R(2,3)*z;
%     w = R(3,1)*x + R(3,2)*y + R(3,3)*z;
%     
%     [lam, th] = cart2sph(u, v, w);
%     
%     y = feval(f, lam, pi/2-th);
% end
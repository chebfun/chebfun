function g = rotate(f, phi, theta, psi)
%ROTATE   Rotates a BALLFUN using Euler angles.
%   Y = ROTATE(F, PHI, THETA, PSI) rotates F using Euler angles phi, theta, 
%   and psi with the ZXZ convention:  Rotate first about the z-axis by an
%   angle phi, then about the (orginal) x-axis by an angle 0<=theta<=pi, 
%   then about new z-axis by an angle psi. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    return
elseif ( nargin == 2 )
    theta = 0;
    psi = 0;
elseif ( nargin == 3 )
    psi = 0;
end

[m, n, p] = size( f );

% f is bandlimited of degree (n,m) so its rotation will be bandlimited will
% limit at most max(n,m). This follows from spherical harmonic theory as
% the rotation of a spherical harmonic of degree l is of degree l, only the
% order will change. The degree l is given by the bandlimit m, while the
% order is given by the bandlimit of n.
p = max(n,p);             % Using this sampling we should exactly recover the rotated f.
n = p + mod(p,2);         % Number of columns must be even.
% Set the number of rows in the sampled grid equal to n/2+2 so that the 
% doubled up grid will have n rows (n/2+2 because the pole is included in
% the sampled grid, and the doubled up grid does not contain -pi).
p = ceil(p/2)+2;
m = ceil(m/2)+2;

r = chebpts(m);

G = zeros(n,p,m);

for k = 1:m
    % Rotate each sphere of radius r(k)
    fsph = extract_spherefun(f, r(k));
    fi = rotate(fsph, phi, theta, psi);
    G(:,:,k) = coeffs2(fi,n,p);
end

G = permute(G,[3,1,2]);
% Convert Cheb values to coeffs
for k = 1:p
    G(:,:,k) = chebtech2.vals2coeffs(G(:,:,k));
end

g = ballfun(G,'coeffs');
end

function g = extract_spherefun(f, r)
% EXTRACT_SPHEREFUN SPHEREFUN corresponding to the value of f at radius r
%   EXTRACT_SPHEREFUN(f, r) is the SPHEREFUN function 
%   g(lambda, theta) = f(r, :, :)

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.coeffs;
[m,n,p] = size(f);

if m == 1
    G = reshape(F(1,:,:),n,p);
else
    % Chebyshev functions evaluated at r
    T = zeros(1,m);
    T(1) = 1; T(2) = r;
    for i = 3:m
        T(i) = 2*r*T(i-1)-T(i-2);
    end

    % Build the array of coefficient of the spherefun function
    G = zeros(n,p);
    for i = 1:p
        G(:,i) = T*F(:,:,i);
    end
end
% Build the spherefun function; coeffs2spherefun takes the theta*lambda matrix
% of coefficients
g = spherefun.coeffs2spherefun(G.');
end

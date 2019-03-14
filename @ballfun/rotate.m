function f = rotate(f, phi, theta, psi)
%ROTATE   Rotates a BALLFUN using Euler angles
%   Y = ROTATE(F, PHI, THETA, PSI) rotates F using Euler angles phi, theta, 
%   and psi with the ZXZ convention:  Rotate first about the z-axis by an
%   angle phi, then about the (orginal) x-axis by an angle 0<=theta<=pi, 
%   then about new z-axis by an angle psi. 
%
% See also spherefun.rotate, diskfun.rotate.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
    return
end

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
P = max(n,p);             % Using this sampling we should exactly recover the rotated f.
N = P + mod(P,2);         % Number of columns must be even.
% Set the number of rows in the sampled grid equal to n/2+2 so that the 
% doubled up grid will have n rows (n/2+2 because the pole is included in
% the sampled grid, and the doubled up grid does not contain -pi).
P = ceil(P/2)+2;
P = P + mod(P,2);

% Sampling grid.
[lam,th] = meshgrid(trigpts(N,[-pi pi]),linspace(0,pi,P));
lam(1,:) = 0;
lam(P,:) = 0;
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

% Get size of evaluation points and Fourier coefficients:
m = m + (1-mod(m,2));
n = n + (1-mod(n,2));
p = p + (1-mod(p,2));
F = coeffs3(f, m, n, p);

% Convert F to chebvals * trigcoeffs
for k = 1:p
   F(:,:,k) = chebtech2.coeffs2vals(F(:,:,k));
end

F = F(ceil(m/2):end,:,:);

% Restrict to [0,1]
m = ceil(m/2);

% Permute F
F = permute(F, [2,3,1]);

% Create empty tensors
G = zeros(N,P,m);

% Loop over r and rotate each spheroid
for k = 1:m
    Fk = F(:,:,k);
    G(:,:,k) = fastSphereEval(Fk.', lam, th).';
end

% Permute G
G = permute(G, [3,1,2]);

% Create ballfun from values
f = ballfun(G);
end

function vals = fastSphereEval(f, lambda, theta)
%FASTSPHEREEVAL     Fast evaluation of a spherefun
%   VALS = FASTSPHEREEVAL(F, LAMBDA, THETA) evaluates the spherefun F at the
%   points (LAMBDA, THETA), where LAMBDA is longitude, with -pi<=LAMBA<=pi, and
%   THETA is co-latitude, with 0<=THETA<=PI. Returns the values of F at these
%   points in the array VALS.
%
%   This is a fast transform analogue to the direct algorithm evaluation of F
%   using Horner's rule. The evaluation costs O(m*n*log(m*n) + N) operations.
%
%   The algorithm in this MATLAB script is based on the paper:
%    [1] D. Ruiz--Antoln and A. Townsend, "A nonuniform fast Fourier transform
%    based on low rank approximation", submitted, 2017.
%
% See also SPHEREFUN/FEVAL and SPHEREFUN/ROTATE.

% Copyright 2017 by The University of Oxford and The CHEBFUN Developers.
% See http://www.chebfun.org/ for CHEBFUN information.

% Primary author: Alex Townsend, January 2017.

% Working accuracy: 
tol = 10*eps;

% Convert to NUFFT2D convention (from Spherefun's convention):
lambda = -lambda/(2*pi);
theta  = -theta/(2*pi);

% Get size of evaluation points and Fourier coefficients:
[M, N] = size(lambda);
[m,n] = size(f);

% Low rank approximation for the coefficients of F:
C = f;
D = eye(n,n);
R = eye(n,n);
rk = size(D, 1);

% Make columns:
lambda = lambda(:);
theta  = theta(:);

% Assign each evaluation point to its closest point on a 2D equispaced grid:
xj = round(n*lambda)/n;
xt = mod(round(n*lambda),n)+1;
yj = round(m*theta)/m;
yt = mod(round(m*theta),m)+1;

% Fourier modes:
mm = -floor(m/2):floor(m/2);
nn = -floor(n/2):floor(n/2);

%%%%%%%%%%%%%%%% FAST TRANSFORM %%%%%%%%%%%%%%%%
% Ay = exp(-2*pi*1i*m*(y(:)-yj(:))*mm/m);
% Find low rank approximation to Ay = U1*V1.':
er = m*(theta-yj);
gam = norm(er, inf);

% Faster version of 
% K1 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
xi = log(log(10/tol)/gam/7);
lw = xi - log(xi) + log(xi)/xi + .5*log(xi)^2/xi^2 - log(xi)/xi^2;
K1 = ceil(5*gam*exp(lw));

% Low rank factors for Ay: 
Dy = diag((1i).^(0:K1-1));
invDy = diag((1./1i).^(0:K1-1));
U1 = real( chebT(K1-1,er/gam) * besselCoeffs(K1, gam) * Dy );
V1 = chebT(K1-1, 2*mm'/m) * invDy;

% Ax = exp(-2*pi*1i*n*(x(:)-xj(:))*nn/n);
% Find low rank approximation to Ax = U2*V2.':
er = n*(lambda-xj);
gam = norm(er, inf);

% Faster version of 
% K2 = ceil(5*gam*exp(lambertw(log(10/tol)/gam/7)));
xi = log(log(10/tol)/gam/7);
lw = xi - log(xi) + log(xi)/xi + .5*log(xi)^2/xi^2 - log(xi)/xi^2;
K2 = ceil(5*gam*exp(lw));

% Low rank factors for Ax:
Dx = diag((1i).^(0:K2-1));
invDx = diag((1./1i).^(0:K2-1));
U2 = real( chebT(K2-1,er/gam) * besselCoeffs(K2, gam) * Dx );
V2 = chebT(K2-1, 2*nn'/n) * invDx;

% Business end of the transform. (Everything above could be considered
% precomputation.)
FFT_cols = zeros(m,rk,K1);
for s = 1:K1
    % Transform in the "col" variable:
    CV1 = bsxfun(@times, C, V1(:,s));
    FFT_cols(:,:,s) = fft( ifftshift(CV1, 1) );
end
FFT_rows = zeros(rk,n,K2);
for r = 1:K2
    % Transform in the "row" variable:
    RV2 = bsxfun(@times, R, V2(:,r)).';
    FFT_rows(:,:,r) = D*fft( ifftshift(RV2,2), [], 2 );
end

% Permute:
XX = permute(FFT_rows, [1 3 2]);
YY = permute(FFT_cols, [3 2 1]);
% Convert to cell for speed:
X = cell(n,1);
for k = 1:n
    X{k} = XX(:,:,k);
end
Y = cell(size(YY,3),1);
for k = 1:m
    Y{k} = YY(:,:,k);
end

% Spread the love out from equispaced points to actual evaluation points. Do
% this K1*K2 times for an accurate transform:
[ii, jj] = ind2sub([m,n], yt(:)+m*(xt(:)-1));
% Do only unique multiplications:
[c, ~, ic] = unique([ii, jj], 'rows');
Y = Y(c(:,1));
X = X(c(:,2));

% Determine the (ii,jj) values that do the same multiplications involving
% Y and X so we can vectorize these operations.
[srt_ic, idic] = sort(ic);
breaks = find(diff(srt_ic)); 
breaks(end+1) = numel(idic); % Include the endpoint
cnt = 1;
ov = ones(K2,1);
vals = zeros(M, N); % Allocate storage for output.
for idx = 1:size(c,1)
    % Get the (ii,jj) values that use the same X and Y.
    kk = idic(cnt:breaks(idx));
    % Do the inner product:
    A = U1(kk,:)*Y{idx};
    vals(kk) = (U2(kk,:).*(A*X{idx}))*ov;  % This is faster than sum2
    cnt = breaks(idx)+1;
end
end

function cfs = besselCoeffs(K, gam)
% The bivarate Chebyshev coefficients for the function f(x,y) = exp(-i*x.*y) on
% the domain [-gam, gam]x[0,2*pi] are given by Lemma A.2 of Townsend's DPhil
% thesis.
arg = -gam*pi/2;
[pp, qq] = meshgrid(0:K-1);
cfs = 4*(1i).^qq.*besselj((pp+qq)/2,arg).*besselj((qq-pp)/2, arg);
cfs(2:2:end,1:2:end) = 0;
cfs(1:2:end,2:2:end) = 0;
cfs(1,:) = cfs(1,:)/2;
cfs(:,1) = cfs(:,1)/2;
end

function T = chebT(n, x)
% Evaluate Chebyshev polynomials of degree 0,...,n at points in x. Use the
% three-term recurrence relation:
N = size(x,1);
T = zeros(N, n+1);
T(:,1) = 1;
if ( n == 1 )
    return
end
T(:,2) = x;
twoX = 2*x;
for k = 2:n
    T(:,k+1) = twoX.*T(:,k) - T(:,k-1);
end
end
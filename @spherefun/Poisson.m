function u = Poisson(f, const, m, n)
%POISSON   Fast Poisson solver for the sphere.
%   POISSON(F, CONST, N) solves
%
%     sin(th)^2u_{th,th} + sin(th)cos(th)u_th + u_{lam,lam} = sin(th)^2*f
%
%   on the unit sphere written in spherical coordinates (lam, th)
%   with integral condition sum2(u) = CONST with a discretization of size N x N.
%
%   POISSON(F, CONST, M, N) is the same as POISSON(F, CONST, N), but with a
%   discretization of size M x N.
%
% EXAMPLE:
%   f = @(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th);
%   exact = @(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
%             sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th);
%   const = 0;
%   u = spherefun.Poisson(f, const, 100);
%   norm(spherefun(exact) - u)

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPERS NOTE:
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in
% the theta- and lambda-direction.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation decouples into n
% linear systems. This form banded matrices.
%
% SOLVE COMPLEXITY:    O(M*N)  with M*N = total degrees of freedom
%
% This is designed to be almost a standalone script to reduce overhead 
% so that it can be used for the numerical simulation of time-dependent 
% PDEs, where this comment is executed hundreds of times. 

if ( nargin < 4 )
    n = m;
end

% Construct useful spectral matrices:
% Please note that DF1m is different than trigspec.diff(m,1) because we 
% take the coefficient space point-of-view and set the (1,1) entry to be 
% nonzero.
DF1m = (1i)*spdiags((-floor(m/2):1:ceil((m-2)/2))', 0, m, m);
DF2m = (1i)^2*spdiags((-floor(m/2):1:ceil((m-2)/2))', 0, m, m).^2;
DF2n = (1i)^2*spdiags((-floor(n/2):1:ceil((n-2)/2))', 0, n, n).^2;

% Multiplication for sin(theta).*cos(theta):
Mcossin = spdiags(.25i*[-ones(m, 1) ones(m, 1)], [-2 2], m, m); 

% Multiplication for sin(theta)^2:
Msin2 = spdiags(.5*[-.5*ones(m, 1) ones(m, 1) -.5*ones(m, 1)], [-2 0 2], m, m);
Im = speye(m);
scl = diag(DF2n); 

% There is a factor of 4 speed up here, by taking account of real 
% solution, and even/odd symmetry.

% Underlying discretization grid:
lam0 = linspace(-pi, pi, n+1)'; 
lam0(end) = [];
th0 = linspace(-pi, pi, m+1)'; 
th0(end) = [];

% Forcing term:
if ( isa(f, 'function_handle') )
    [rhs_lam, rhs_theta] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_theta);
    F = trigtech.vals2coeffs(F);
    F = Msin2*trigtech.vals2coeffs(F.').';
elseif ( isa(f, 'spherefun') )
    F = Msin2*coeffs2(f, n, m);
elseif ( isa( f, 'double' ) )
    F = Msin2*f;       % Get trigcoeffs2 of rhs.
end

% Matrix for solution's coefficients:
CFS = zeros(m, n);

% Form discretization of the theta-dependent operator:
L = Msin2*DF2m + Mcossin*DF1m;

% Solve for the even modes:
k_even = [floor(n/2)-1:-2:1 floor(n/2)+3:2:n];
for k = k_even
    CFS(:,k) = (L + scl(k)*Im) \ F(:,k);
end

% Solve for the odd modes:
k_odd = [floor(n/2):-2:1 floor(n/2)+2:2:n];
for k = k_odd
    CFS(:,k) = (L + scl(k)*Im) \ F(:,k);
end

% Now do the equation where we need the integral constraint:
% We will take X_{n/2+1,:} en = const.

% First, let's project the rhs to have mean zero:
k = floor(n/2) + 1;
floorm = floor(m/2);
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
ii = [1:floorm floorm+2:m];
F(floor(m/2)+1,k) = -en(ii)*F(ii, k)./en(floor(m/2)+1);

% Second, solve: 
CFS(:, k) = [ en ; L( ii, :) ] \ [ const ; F(ii, k) ];
u = spherefun.coeffs2spherefun(CFS); 

end
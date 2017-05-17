function u = poisson(f, const, m, n)
%POISSON   Fast Poisson solver for the sphere.
%   POISSON(F, C, N) solves laplacian(U) = F on the unit sphere, which in 
%   spherical coordinates (lam, th) is
%
%     sin(th)^2U_{th,th} + sin(th)cos(th)U_th + U_{lam,lam} = sin(th)^2*F
%
%   The equation is discretized on an N x N grid in spherical coordinates.
%   The integral of F is assumed to be zero, which is the compatibility
%   constraint for there to exist a solution to the Poisson problem on the
%   sphere. The mean value of the solution U is set to C.  This
%   function returns a SPHEREFUN representing the solution.
%
%   POISSON(F, C, M, N) is the same as POISSON(F, C, N), but with a
%   discretization of size M x N.
%
% EXAMPLE:
%   f = @(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th);
%   exact = @(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
%             sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th);
%   u = spherefun.poisson(f, 0, 100);
%   norm(spherefun(exact) - u)
%   mean2(u)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
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

% TODO: 
% Make this code adaptive! 

if ( nargin < 4 )
    n = m;
end

% Construct useful spectral matrices:
% Please note that DF1m is different than trigspec.diff(m,1) because we 
% take the coefficient space point-of-view and set the (1,1) entry to be 
% nonzero.
DF1m = trigspec.diffmat(m, 1, 1);
DF2m = trigspec.diffmat(m, 2); 
DF2n = trigspec.diffmat(n, 2); 

% Multiplication for sin(theta).*cos(theta):
% Below is equivalent to 
% Mcossin = spdiags(.25i*[-ones(m, 1) ones(m, 1)], [-2 2], m, m); 
cfs = trigtech(@(theta) sin(pi*theta).*cos(pi*theta));
Mcossin = trigspec.multmat(m, cfs.coeffs); 

% Multiplication for sin(theta)^2:
% Below is equivalent to
% Msin2 = spdiags(.5*[-.5*ones(m, 1) ones(m, 1) -.5*ones(m, 1)], [-2 0 2], m, m);
cfs = trigtech(@(theta) sin(pi*theta).^2);
Msin2 = trigspec.multmat(m, cfs.coeffs);
Im = speye(m);
scl = diag(DF2n); 

% There is a factor of 4 speed up here, by taking account of real 
% solution, and even/odd symmetry.

% Underlying discretization grid:
lam0 = trigpts(n,[-pi, pi]); 
th0 = trigpts(m,[-pi, pi]); 

% Forcing term:
if ( isa(f, 'function_handle') )
    [rhs_lam, rhs_theta] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_theta);
    tol = 1e5*max(abs(F(:)))*chebfunpref().cheb2Prefs.chebfun2eps;
    F = trigtech.vals2coeffs(F);
    F = trigtech.vals2coeffs(F.').';
elseif ( isa(f, 'spherefun') )
    tol = 1e5*vscale(f)*chebfunpref().cheb2Prefs.chebfun2eps;
    F = coeffs2(f, n, m);
elseif ( isa( f, 'double' ) )
    tol = 1e5*chebfunpref().cheb2Prefs.chebfun2eps;
    F = f;       % Get trigcoeffs2 of rhs.
end

% First, let's project the rhs to have mean zero:
k = floor(n/2) + 1;
floorm = floor(m/2);
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
ii = 1:m;
meanF = en(ii)*F(ii, k)/en(floor(m/2)+1);

% Check that the mean of F is zero (or close enough).  If it is not then
% issue a warning
if ( abs(meanF) > tol )
    warning('CHEBFUN:SPHEREFUN:POISSON:meanRHS',...
       ['The integral of the right hand side may not be zero, which is '...
        'required for there to exist a solution to the Poisson '...
        'equation. Subtracting the mean off the right hand side now.']);
end        
F(floor(m/2)+1,k) = F(floor(m/2)+1,k)-meanF;

% Multiply the right hand side by (sin(theta)).^2
F = Msin2*F;

% Matrix for solution's coefficients:
CFS = zeros(m, n);

% Form discretization of the theta-dependent operator:
L = Msin2*DF2m + Mcossin*DF1m;

% Solve for the even modes:
kk = [floor(n/2):-1:1 floor(n/2)+2:n];
%k_even = [floor(n/2)-1:-2:1 floor(n/2)+3:2:n];
for k = kk
    CFS(:,k) = (L + scl(k)*Im) \ F(:,k);
end

% % Solve for the odd modes:
% k_odd = [floor(n/2):-2:1 floor(n/2)+2:2:n];
% for k = k_odd
%     CFS(:,k) = (L + scl(k)*Im) \ F(:,k);
% end

% Now do the equation where we need the integral constraint:
% We will take X_{n/2+1,:} en = 0.

% Second, solve: 
k = floor(n/2) + 1;
ii = [1:floorm floorm+2:m];
CFS(:, k) = [ en ; L( ii, :) ] \ [ 0 ; F(ii, k) ];
u = spherefun.coeffs2spherefun( CFS ) + const; 

end
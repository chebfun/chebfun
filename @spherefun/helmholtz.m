function u = helmholtz(f, K, m, n)
%HELMHOLTZ   Fast Helmholtz solver for the sphere.
%    U = HELMHOLTZ(F, K, N) solves U_xx + U_yy + U_zz + K^2U = F on the sphere
%    for U with a discretization of size N x N. F should be a SPHEREFUN and the
%    solution is returned as a SPHEREFUN.
%
%    HELMHOLTZ(F, K, M, N) same as HELMHOLTZ(F, K, N), but with a
%    discretization of size M x N.
%
%  Example:
%    K = 100; m = 1000; n = m;
%    f = spherefun( @(x,y,z) cos(x.*y.*z) );
%    u = spherefun.helmholtz(f, K, m, n);
%    plot( u )

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPERS NOTE:
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in
% the theta- and lambda-direction.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation decouples into N
% linear systems. This form banded matrices.
%
% SOLVE COMPLEXITY:    O(M*N)  with M*N = total degrees of freedom

% Solve standard Helmholtz equation. This parameter is kept for developers.
c = 1;

% If the call is helmholtz(f, K, m), then set n
if ( nargin < 4 )
    n = m;
end

if ( K == 0 )
    u = spherefun.poisson(f, 0, m, n);
    return
end

% Check for eigenvalues:
e = eig( [ -1 K^2 ; 1 0 ] ); % Is K = sqrt(l*(l+1)), where l is an integer?
e = e(e>0); e = e(abs( e - round(e) ) < 1e-13 );
if ( ~isempty( e ) )
    error('SPHEREFUN:HELMHOLTZ:EIGENVALUE',...
            'There are infinitely many solutions since K is an eigenvalue of the Helmholtz operator.')
end

% If m or n are non-positive then throw an error
if ( ( m <= 0 ) || ( n <= 0 ) )
    error('CHEBFUN:SPHEREFUN:HELMHOLTZ:badInput',...
        'Discretization sizes should be positve numbers');
end

% If m and n are 1, the solution is easy.  It's just a constant given by
% 1/K^2*mean(F).
if ( ( m == 1 ) && (n == 1 ) )
    f = spherefun(f);
    u = mean2(f)/K^2 + 0*f;
    return;
end

% Make m even so that the pole at theta=0 is always sampled.
m = m + mod(m,2);

% Construct useful spectral matrices:
Im = speye(m);

% Please note that DF1m here is different than trigspec.diff(m,1) because we
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

% Underlying discretization grid:
lam0 = linspace(-pi, pi, n+1)';
lam0( end )=[];
th0 = linspace(-pi, pi, m+1)';
th0( end )=[];

% Forcing term:
if ( isa(f, 'function_handle') )
    [rhs_lam, rhs_th] = meshgrid(lam0, th0);
    F = feval(f, rhs_lam, rhs_th);      % Get (trigvals,trigvals) of rhs
    F = trigtech.vals2coeffs(F);        % Get in Fourier basis
    F = trigtech.vals2coeffs(F.').';
elseif ( isa(f, 'spherefun') )
    F = coeffs2(f, n, m);
end

% Calculate the integral constraint constant:
k = floor(n/2)+1;
floorm = floor(m/2);
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
int_const = en*F(:,k)/K^2;

% Multiple rhs by sin(th)^2 and divide by K^2:
F = Msin2 * F / K^2;

% Want to solve
%    X L^T + X DF^T = F
% subject to zero integral constraints, i.e.,  w^T X_0  = 0.

% Note that the matrix equation decouples because DF is diagonal.

% Solve decoupled matrix equation for X, one row at a time:
CFS = zeros(m, n);
L = c*(Msin2*DF2m + Mcossin*DF1m)/K^2 + Msin2;
scl = c*diag(DF2n)/K^2;
for k = [floor(n/2):-1:1 floor(n/2)+2:n]
    CFS(:,k) = (L + scl(k)*Im) \ F(:,k);
end

% Now, do zeroth mode:
k = floor(n/2)+1;
ii = [1:floorm floorm+2:m];
CFS(:, k) = [ en ; L( ii, :) ] \ [ int_const ; F(ii, k) ];

% Now, convert to a spherefun object:
u = spherefun.coeffs2spherefun(CFS);

end
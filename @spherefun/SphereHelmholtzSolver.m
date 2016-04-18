function CFS = Helmholtz( f, K, m, n )
% HELMHOLTZ  Fast Helmholtz solver for the sphere
%
% HELMHOLTZ( f, K, n ) solves u_xx+u_yy+u_zz +K^2u = f on the sphere with
% a discretization of size nxn.
%
% HELMHOLTZ( f, K, m, n) samee as HELMHOLTZ( f, K, n ), but with a
% discretization of size m x n.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPERS NOTE:
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in
% the theta- and lambda-direction.
%
% LINEAR ALGEBRA: Matrix equations.  The matrix equation decouples into n
% linear systems. This form banded matrices.
%
% SOLVE COMPLEXITY:    O( m*n )  N = m*n = total degrees of freedom

% Solve standard Helmholtz equation.  This parameter is kept for developers.
c = 1;  

if ( K == 0 )
    CFS = SpherePoissonSolver(f, m, n, K);
    return
end

% Construct useful spectral matrices:
Im = speye(m);
DF1m = trigspec.diffmat( m, 1);    % 1st order Fourier diffmat
DF2m = trigspec.diffmat( m, 2);  % 2nd order Fourier diffmat
DF2n = trigspec.diffmat( n, 2);
Mcossin = spdiags( .25i*[-ones(m,1) ones(m,1)], [-2 2], m ,m ); % multiplication for sin(theta).*cos(theta)
Msin2 = spdiags( .5*[-.5*ones(m,1) ones(m,1) -.5*ones(m,1)], [-2 0 2], m ,m );% multiplication for sin(theta)^2

% Discretization sizes:
lam0 = pi*trigpts( n );     %
th0 = pi*trigpts( m );      % GRID

% Forcing term:
if ( isa(f, 'function_handle') || isa(f, 'spherefun') )
    [rhs_lam, rhs_th] = meshgrid( lam0, th0 );
    F = feval( f, rhs_lam, rhs_th );      % Get (trigvals,trigvals) of rhs
    F = trigtech.vals2coeffs( F );         % Get in Fourier basis
    F = trigtech.vals2coeffs( F.' ).';
    F = Msin2*F;
elseif ( isa( f, 'double' ) )
    F = Msin2*f;                                 % Get in Fourier basis
end

% Want to solve        X L^T   +    X DF^T = F
%
%   subject to     zero integral constraints, i.e.,  w^T X_0  = 0

% Note that the matrix equation decouples because DF is diagonal.

% Solve decoupled matrix equation for X, one row at a time:
CFS = zeros(m, n);

L = c*(Msin2*DF2m + Mcossin*DF1m)/K^2 + Msin2;

scl = c*diag( DF2n )/K^2;

for k = [floor(n/2):-1:1 floor(n/2)+2:n]
    
    CFS(:,k) = ( L + scl(k)*Im ) \ F(:,k);
    
end

% Now, do zeroth mode:
k = floor(n/2)+1;
floorm = floor(m/2);
mm = (-floorm:ceil(m/2)-1);
en = 2*pi*(1+exp(1i*pi*mm))./(1-mm.^2);
en([floorm, floorm + 2]) = 0;
ii = [1:floorm floorm+2:m];
int_const = en*F(:,k)/K^2;
CFS(:, k) = [ en ; L( ii, :) ] \ [ int_const ; F(ii, k) ];

end
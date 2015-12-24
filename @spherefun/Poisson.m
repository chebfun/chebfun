function u = Poisson( f, int_const, m, n )
% POISSON  Fast Poisson solver for the sphere
% 
% POISSON( f, const, n ) solves  
%  
%  sin(th)^2u_thth   + sin(th)cos(th)u_th +  u_{lam,lam}  =  sin(th)^2 * f   
%
% on the unit sphere written in spherical coordinates (lam, th) 
% with integral condition  sum2( u ) = const with a discretization of 
% size N x N.
%
% POISSON( f, const, m, n ) same as POISSON( f, const, n), but with a
% discretization of size m x n. 
%
% EXAMPLE: 
%  f = @(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th);
%  exact = @(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
%            sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th);
%  int_const = 0;
%  u = spherefun.Poisson( f, int_const, 100);
%  norm( spherefun(exact) - u )

% DEVELOPERS NOTE:
% 
% METHOD: Spectral method (in coeff space). We use the Fourier basis in 
% the theta- and lambda-direction.
%
% LINEAR ALGEBRA: Matrix equations.  The matrix equation decouples into n 
% linear systems. This form banded matrices. 
%
% SOLVE COMPLEXITY:    O( m*n )  N = m*n = total degrees of freedom
%
%  Alex Townsend, July 2015. 

if ( nargin < 4 ) 
    n = m; 
end

% Double up to do DFS method: 
n = 2*n; 

% Discretization grid: 
lam0 = pi*trigpts( m );  
th0 = pi*trigpts( n );      

% Calculate trigcoeffs of the forcing term: 
[ll, tt] = meshgrid( lam0, th0 ); 
vals = sin(tt).^2 .* feval(f, ll, tt );
G = trigtech.vals2coeffs( trigtech.vals2coeffs( vals ).' ).'; 

% Construct useful spectral discretization matrices: 
DF1n = 1i*spdiags((-n/2:n/2-1)', 0, n, n);  % 1st order Fourier diffmat
DF2n = spdiags( -(-n/2:n/2-1)'.^2, 0, n, n);% 2nd order Fourier diffmat
DF2m = spdiags( -(-m/2:m/2-1)'.^2, 0, m, m);
% Multiplication operator for sin(theta):
Msin = spdiags( .25i*[-ones(n,1) ones(n,1)], [-2 2], n ,n ); 
% Multiplication operator for sin(theta).*cos(theta):
Msin2 = spdiags( .5*[-.5*ones(n,1) ones(n,1) -.5*ones(n,1)], [-2 0 2], n ,n );

% Form discretization of the theta-dependent operator: 
L = Msin2*DF2n + Msin*DF1n;

% We want to solve 
%     L * X  +  X * DF2^T  = (sin(th))^2*g  
% such that X(m/2+1, m/2+1) = int_const.  

% Notice that the matrix equation decouples and we can solve for X one
% column at a time: 
CFS = zeros( n, m ); 
I = speye( n );
% We leave out the (m/2+1)th column because that needs a constraint
for k = [1:m/2 m/2+2:m]  
    CFS(:, k) = ( L +  DF2m(k, k)*I ) \ G(:, k);
end

% Now do the equation where we need the integral constraint. 
% We ensure that X_{n/2+1,:}*en = int_const. 
en = 2*pi*(1+exp(1i*pi*(-n/2:n/2-1)))./(1-(-n/2:n/2-1).^2);
en([n/2, n/2 + 2]) = 0;
CFS(:, m/2+1) = [ en ;  L([1:n/2 n/2+2:n], :)] \ ... 
                                 [ int_const ; G([1:n/2 n/2+2:n], m/2+1) ];

% Convert to VALUES on the underlying grid: 
VALS = trigtech.coeffs2vals( trigtech.coeffs2vals( CFS ).' ).'; 

% Now restrict down to region of interest:
VALS = VALS([n/2+1:n 1], :);

% Finally, make a spherefun obxject out of the values: 
u = spherefun( real( VALS ) ); 

% - DEBUG
% % Plot: 
% subplot(1,3,1)
% surf(ll, tt, real(VALS),'edgealpha',0,'facecolor','interp'), 
% ylim([-pi pi])
% norm( imag( VALS ) )
% 
% % Compare to exact solution:
% subplot(1,3,2)
% EXACT = exact(ll,tt);
% EXACT_CFS = trigtech.vals2coeffs( trigtech.vals2coeffs( EXACT ).' ).'; 
% surf(ll, tt ,EXACT ,'edgealpha',0,'facecolor','interp'), 
% ylim([-pi pi])
% 
% subplot(1,3,3)
% surf(ll, tt, abs( VALS - EXACT ),'edgealpha',0,'facecolor','interp' )

end
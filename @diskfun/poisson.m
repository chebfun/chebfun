function u = poisson( f, bc, m, n )
% POISSON    Fast Poisson solver with Dirichlet boundary conditions for the disk. 
%   U = POISSON(F, BC, N) solves laplacian(U) = F on the unit disk, which
%   in polar coordinates (theta, r) is 
%  
%       r^2 u_rr   +   r u_r   +  u_{theta theta}   =   r.^2 * f   
%
%   The equation is discretized on an N x N grid in polar coordinates. The
%   solution is imposed to satisfy the Dirichlet data 
%   U(1,theta) = @(theta) BC(theta). The solution is returned as a 
%   diskfun object. (N must be even.)
%
%   U = POISSON(F, BC, M, N) is the same as POISSON(F, BC, N) but uses a
%   discretization size of M x N. (N must be even.)
%
%   F may be a function handle in polar coordinates, a diskfun, or a matrix 
%   of Fourier-Chebyshev coefficients. 
% 
%   EXAMPLE: 
%    bc = @(th) 0*th;              
%    f = @(th, r) -1 + 0*th;            
%    u = diskfun.poisson( f, bc, 100); 
%
% See also SPHEREFUN/POISSON

% DEVELOPER'S NOTE: 
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in the
% periodic theta-direction and the Ultraspherical spectral method in the
% radial direction.   
%
% LINEAR ALGEBRA: The matrix equation decouples into n linear systems. 
% Using parity properties, these can be expressed as tridiagonal + rank-1 matrices
% and are solved using the Sherman-Morrison (Woodbury) formula. 
%
% SOLVE COMPLEXITY:    O( m*n )  N = m*n = total degrees of freedom
%
%   Alex Townsend, July 2015. (Updated Feb 2016,  Heather Wilber)
%
% For more details, see 
% H. Wilber, Numerical computing with functions on the sphere and disk,
% Master?s thesis, Boise State University, 2016.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if( nargin < 4 )
    n = m;
end

% Call Helmholtz code with K = 0: 
u = diskfun.helmholtz(f, 0, bc, m, n );

end


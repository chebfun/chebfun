function u = poisson( f, g, m, n )
%POISSON   Fast Poisson solver for the rectangle.
%   POISSON(F, G, N) solves laplacian(U) = F on the domain of f with
%   Dirichlet boundary conditions given by G. That is, U satisfies
%
%     U_{x,x} + U_{y,y} = F, on [a,b]x[c,d]    U = G on boundary
%
%   The equation is solved using an N x N discretization. 
%
%   POISSON(F, G, M, N) is the same as POISSON(F, G, N), but with an M x N
%   discretization. 
%
% EXAMPLE:
%   f = chebfun2( @(x,y) 1 + 0*x, [-1 2 0 1]);
%   u = chebfun2.poisson(f, 0, 100);
%   plot(u)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPERS NOTE:
%
% METHOD: Spectral method (in coeff space). We use a C^{(3/2)} basis to
% discretize the equation, resulting in a discretization of the form AX+XA
% = F, where A is a symmetric tridiagonal matrix. 
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation is solved by the
% alternating direction implicit (ADI) method. 
%
% SOLVE COMPLEXITY:  O(M*N*log(MAX(M,N))log(1/eps))  with M*N = total 
% degrees of freedom. 

if ( nargin == 3 ) 
    % Call is POISSON(F, G, N) so employ an NxN discretization:
    n = m;
end

% Compute the C^(3/2) coefficients of f(x,y):
F = Chebyshev2ultra( coeffs2( f, m, n ) );

% Construct M, the multiplication matrix for (1-x^2) in the C^(3/2) basis
j = (0:n-1)';
dsub = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j+2);
dsup = -1./(2*(j+3/2)).*(j+1).*(j+2)*1/2./(1/2+j);
d = -dsub - dsup;
M = spdiags([dsub d dsup], [-2 0 2], n, n);

% Construct D^{-1}, which undoes the scaling from the Laplacian identity
D_inv = spdiags(-1./(j.*(j+3)+2), 0, n, n);

% Construct T = D^{-1} * M and extract T_even and T_odd
T = D_inv * M;
Te = T(1:2:n, 1:2:n);
To = T(2:2:n, 2:2:n);

% Solve the Sylvester equation for the even and odd coefficients
% separately. The goal is to replace lyap with something faster.
rhs = -D_inv * F * D_inv;
Yee = lyap( Te, rhs(1:2:n, 1:2:n) );
Yoo = lyap( To, rhs(2:2:n, 2:2:n) );
Yeo = lyap( Te, To.', rhs(1:2:n, 2:2:n) );
Yoe = lyap( To, Te.', rhs(2:2:n, 1:2:n) );
Y = zeros(n,n);
Y(1:2:n, 1:2:n) = Yee;
Y(2:2:n, 2:2:n) = Yoo;
Y(1:2:n, 2:2:n) = Yeo;
Y(2:2:n, 1:2:n) = Yoe;

% Convert back to Chebyshev
U = ultra1mx2Chebyshev( Y );
u = chebfun2( U, 'coeffs' );

end
function u = poisson( f, g, m, n )
%POISSON   Fast Poisson solver for the rectangle.
%   POISSON(F, G, N) solves laplacian(U) = F on the domain of f with
%   Dirichlet boundary conditions given by G. That is, U satisfies
%
%     U_{x,x} + U_{y,y} = F, on [a,b]x[c,d]    U = G on boundary
%
%   The equation is solved using an N x N discretization. G can be a
%   scalar or any chebfun2 object satisfying the Dirichlet data.
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
% 
% AUTHORS: Dan Fortunato (dan.fortunato@gmail.com ) and Alex 
%          Townsend (townsend@cornell.edu). 

% Solve for u on the same domain as f, adjust diffmat to including scaling:
dom = f.domain;
scl_x = (2/(dom(2)-dom(1)))^2;
scl_y = (2/(dom(4)-dom(3)))^2;

if ( nargin == 3 )
    % Call is POISSON(F, G, N) so employ an NxN discretization:
    n = m;
end

% Compute the Chebyshev coefficients of f(x,y):
F = coeffs2( f, m, n );

% Solver only deals with zero homogeneous Dirichlet conditions. Therefore,
% if nonzero Dirichlet conditions are given, we solve lap(u) = f with u|bc = g
% as    u = v + w, where v|bc = g, and lap(w) = f - lap(v), w|bc = 0:
if ( isa(g, 'double') )
    BC = zeros(m, n);
    BC(1,1) = g;
elseif ( isa(g, 'chebfun2') )
    if ( g.domain ~= f.domain )
        error('CHEBFUN2:POISSON:BC', ...
            'Dirichlet data should be on the same domain as F.');
    else
        BC = coeffs2(g, m, n);
        % Adjust the rhs:
        F = F - coeffs2(lap(g), m, n);
    end
elseif ( isa(g, 'function_handle') )
    g = chebfun2(g, f.domain);
    BC = coeffs2(g, m, n);
    % Adjust the rhs:
    F = F - coeffs2(lap(g), m, n);
else
    error('CHEBFUN2:POISSON',...
        'Dirichlet data needs to be given as a scalar or function')
end

% Convert rhs to C^{(3/2)} coefficients: 
F = Chebyshev2ultra( F );

% Construct M, the multiplication matrix for (1-x^2) in the C^(3/2) basis
jj = (0:n-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mn = spdiags([dsub d dsup], [-2 0 2], n, n);
% Construct D^{-1}, which undoes the scaling from the Laplacian identity
invDn = spdiags(-1./(jj.*(jj+3)+2), 0, n, n);
Tn = scl_y * invDn * Mn;

jj = (0:m-1)';
dsub = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj+2);
dsup = -1./(2*(jj+3/2)).*(jj+1).*(jj+2)*1/2./(1/2+jj);
d = -dsub - dsup;
Mm = spdiags([dsub d dsup], [-2 0 2], m, m);
invDm = spdiags(-1./(jj.*(jj+3)+2), 0, m, m);

% Construct T = D^{-1} * M:
Tm = scl_x * invDm * Mm;
F = invDm * F * invDn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Alternating Direction Implicit method %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve TmX + XTn' = F using ADI, which requires O(n^2log(n)log(1/eps)) 
% operations:

% Calculate ADI shifts based on bounds on the eigenvalues of Tn and Tm:
a = -4/pi^2 * scl_y;
b = -39*n^-4 * scl_y;
c = 39*m^-4 * scl_x;
d = 4/pi^2 * scl_x;
[p, q] = ADIshifts(a, b, c, d, 1e-14);

% Run the ADI method:
X = zeros(m, n);
A = Tm; B = -Tn';
Im = speye(m);
In = speye(n);
for j = 1:numel(p)
    X = (F-(A+q(j)*Im)*X) / (B+q(j)*In);
    X = (A+p(j)*Im) \ ( F - X*(B+p(j)*In) );
end

% Convert back to Chebyshev
X = ultra1mx2Chebyshev( X );
X = X + BC;
u = chebfun2( X, f.domain, 'coeffs' );

end

function [p, q] = ADIshifts(a, b, c, d, tol)
% ADISHIFTS  ADI shifts for AX-XB=F when the eigenvalues of A (B) are in [a,b] and
% the eigenvalues of B (A) are in [c,d]. WLOG, we require that a<b<c<d and 0<tol<1.
gam = (c-a)*(d-b)/(c-b)/(d-a);                 % Cross-ratio of a,b,c,d
% Calculate Mobius transform T:{-alp,-1,1,alp}->{a,b,c,d} for some alp:
alp = -1 + 2*gam + 2*sqrt(gam^2-gam);          % Mobius exists with this t
A = det([-a*alp a 1; -b b 1 ; c c 1]);         % Determinant formulae for Mobius
B = det([-a*alp -alp a; -b -1 b ; c 1 c]);
C = det([-alp a 1; -1 b 1 ; 1 c 1]);
D = det([-a*alp -alp 1; -b -1 1; c 1 1]);
T = @(z) (A*z+B)./(C*z+D);                     % Mobius transfom
J = ceil( log(16*gam)*log(4/tol)/pi^2 ); % No. of ADI iterations
if ( alp > 1e7 )
    K = (2*log(2)+log(alp)) + (-1+2*log(2)+log(alp))/alp^2/4;
    m1 = 1/alp^2; 
    u = (1/2:J-1/2)*K/J; 
    dn = sech(u) + .25*m1*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u); 
else
    K = ellipke( 1-1/alp^2 );
    [~, ~, dn] = ellipj((1/2:J-1/2)*K/J,1-1/alp^2);% ADI shifts for [-1,-1/t]&[1/t,1]
end
p = T( -alp*dn ); q = T( alp*dn );             % ADI shifts for [a,b]&[c,d]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERSION CODES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = Chebyshev2ultra( X )
% Convert a matrix of Chebyshev coefficients to a matrix of C^(3/2)
% coefficients.

% First convert the matrix of Chebyshev coefficients to a matrix of
% Legendre coefficients:
[m, n] = size(X); 
if ( m == n && n <= 4000 )
    S = cheb2leg_mat( n );
    X = (S * X) * S.';
elseif ( max( m, n) <= 4000 )
    Sn = cheb2leg_mat( n );
    Sm = cheb2leg_mat( m );
    X = (Sm * X) * Sn.'; 
else
    X = cheb2leg( cheb2leg( X ).' ).';
end

% Now, convert the matrix of Legendre coefficient to a matrix of
% ultraspherical coefficients:
S1 = Legendre2ultraMat( size(X,1) );
S2 = Legendre2ultraMat( size(X,2) );
X = S1 * X * S2.';

end

function X = ultra1mx2Chebyshev( X )
% Convert a matrix of (1-x^2)(1-y^2)C^(3/2)(x)C^(3/2)(y) coefficients
% to a matrix of Chebyshev coefficients.

% First, convert the matrix of (1-x^2)(1-y^2)C^(3/2)(x)C^(3/2)(y) coefficients
% to Legendre coefficients:
[m, n] = size( X );
S1 = ultra1mx2Legendre( m );
S2 = ultra1mx2Legendre( n );
X = S1 * X * S2.';

% Now, convert the matrix of Legendre coefficient to a matrix of Chebyshev
% coefficients:
if ( m == n && n <= 4000 )
    S = leg2cheb_mat( n );
    X = (S * X) * S.';
elseif ( max( m, n) <= 4000 )
    Sn = leg2cheb_mat( n );
    Sm = leg2cheb_mat( m );
    X = (Sm * X) * Sn.'; 
else
    X = leg2cheb( leg2cheb( X ).' ).';
end
end

function S = Legendre2ultraMat( n )
% Conversion matrix from Legendre coefficients to C^(3/2).

lam = 1/2;
dg = lam./(lam + (2:n-1))';
v = [1 ; lam./(lam+1) ; dg];
w = [0 ; 0 ; -dg];
S = spdiags( [v w], [0 2], n, n );

end

function S = ultra1mx2Legendre( n )
% Conversion matrix from (1-x^2)C^(3/2) to Legendre.

d = ones(n, 1);
S = spdiags(((1:n).*(2:(n+1))./2./(3/2:n+1/2))', 0, n, n);
S = spdiags( [d,-d], [0,-2], n, n ) * S;

end

function L = cheb2leg_mat( N ) 
% Construct the cheb2leg conversion matrix.

% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

L = zeros(N, N); 
for j = 0:N-1
    for k = j+2:2:N-1
        L(j+1, k+1) = -k*(j+.5)*(vals((k-j-2)+1)./(k-j)).*(vals((k+j-1)+1)./(j+k+1));
    end
end
c = sqrt(pi)/2;
for j = 1:N-1
    L(j+1, j+1) = c./vals( 2*j+1 ); 
end
L(1,1) = 1; 

end 

function M = leg2cheb_mat( N )
% Construct the leg2cheb conversion matrix.

% This for-loop is a faster and more accurate way of doing:
% Lambda = @(z) exp(gammaln(z+1/2) - gammaln(z+1));
% vals = Lambda( (0:2*N-1)'/2 );
vals = zeros(2*N,1);
vals(1) = sqrt(pi);
vals(2) = 2/vals(1);
for i = 2:2:2*(N-1)
    vals(i+1) = vals(i-1)*(1-1/i);
    vals(i+2) = vals(i)*(1-1/(i+1));
end

M = zeros(N, N); 
for j = 0:N-1
    for k = j:2:N-1
        M(j+1, k+1) = 2/pi*vals((k-j)+1).*vals((k+j)+1);
    end
end
M(1,:) = .5*M(1,:); 

end
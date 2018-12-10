function u = poisson_lr( f,use_fiadi, varargin)
%POISSON   Fast Poisson solver for the rectangle.
%   POISSON(F) solves laplacian(U) = F on the domain of F with zero
%   Dirichlet boundary conditions. That is, U satisfies
%
%     U_{x,x} + U_{y,y} = F, on [a,b]x[c,d], with U = 0 on boundary
%
%   The equation is solved using an adaptively determined discretization
%   size.
%
%   POISSON(F, G) solves using Dirichlet boundary conditions given by G. G
%   can be a scalar, a function handle, or any chebfun2 object satisfying
%   the Dirichlet data.
%
%   POISSON(F, G, N) is the same as POISSON(F, G), but uses an N x N tensor
%   product discretization to solve the equation.
%
%   POISSON(F, G, M, N) is the same as POISSON(F, G, N), but with an M x N
%   tensor product discretization, where N is the number of coeffcieints in
%   the x-direction and M is the number in the y-direction.
%
% EXAMPLE:
%   f = chebfun2( @(x,y) 1 + 0*x, [-1 2 0 1]);
%   u = chebfun2.poisson(f);
%   plot(u)
%
% See also DISKFUN/POISSON, SPHEREFUN/POISSON.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPER'S NOTE:
%
% METHOD: Spectral method (in coefficient space). We use a C^{(3/2)} basis
% to discretize the equation, resulting in a discretization of the form
% AX + XA = F, where A is a symmetric tridiagonal matrix.
%
% LINEAR ALGEBRA: Matrix equations. The matrix equation is solved in 
% low rank form by the factored alternating direction implicit (fADI) 
% method or factored-independent (FI-ADI) method. 
%
% SOLVE COMPLEXITY:  O(M*N*log(MAX(M,N))*log(1/eps)) with M*N = total
% degrees of freedom.
%
% AUTHORS: Dan Fortunato (dan.fortunato@gmail.com)
%          Alex Townsend (townsend@cornell.edu)
%          Heather Wilber (heatherw3521@gmail.com)
%
% The fast Poisson solver is based on:
%
% D. Fortunato and A. Townsend, Fast Poisson solvers for spectral methods,
% Submitted, 2017.
%
% Add fi-adi paper. 

% Solve for u on the same domain as f, adjust diffmat to include scaling:
dom = f.domain;
scl_x = (2/(dom(2)-dom(1)))^2;
scl_y = (2/(dom(4)-dom(3)))^2;

if ( nargin == 1 )
    % Call is POISSON(F) so set G = 0 and use chebop2 to determine the
    % discretization size:
    g = 0;
    N = chebop2(@(u) diff(u,2,1) + diff(u,2,2), f.domain);
    N.bc = g;
    u = N \ f;
elseif ( nargin == 2 )
    % Call is POISSON(F, G) so use chebop2 to determine the discretization
    % size:
    g = varargin{1};
    if ( ~isa(g, 'chebfun2') )
        g = chebfun2(g, f.domain);
    end
    N = chebop2(@(u) diff(u,2,1) + diff(u,2,2), f.domain);
    N.lbc = feval(g, f.domain(1), ':');
    N.rbc = feval(g, f.domain(2), ':');
    N.dbc = feval(g, ':', f.domain(3));
    N.ubc = feval(g, ':', f.domain(4));
    u = N \ f;
else
    % We are given a discretization size, so no need to use chebop2
    %g = varargin{1};
    g = 0; %for now zero bc.
    m = varargin{1};
    n = varargin{2};
    %if ( nargin == 3 )
        % Call is POISSON(F, G, N) so employ an NxN discretization:
       % n = m;
    %else
        % Call is POISSON(F, G, M, N):
    %    n = varargin{3};
    %end

    % Set the error tolerance for ADI
    tol = chebfun2eps();

   %% TO DO CODE:
    %TODO: add an adaptive check for fADI or fi-ADI: 
    %determine whether we should use fADI or fi-ADI: 
    % this is a simple strategy that checks whether decay is 
    % approximately geometric. If it is and the rank is high enough, 
    % we use fi-ADI. Otherwise, we just use fADI.  
        %pivots should be positive and ordered large to small
        
    %for now, we will use flags because we need to test the 
    % crossover point in the rank to make a rule of thumb. 
 %%  
    % Set up RHS for low rank ADI:  
    [U, S, V] = cdr(f);  
     U = chebcoeffs(U, m); 
     V = chebcoeffs(V, n); 
     S = diag(S); 
     sgn = sign(S); 
     S = abs(S);
     
   if use_fiadi ==1 
       [S, ids] = sort(S, 'descend');
       U = U(:, ids)*diag(sgn(ids)); 
       S = diag(S); 
       V = V(:, ids); 
   else
       U = U*diag(sgn)*diag(sqrt(S)); 
       V = V*sqrt(diag(S)); 
   end
     
  %% TO-DO: low rank RHS non-zero BC adjustments
  %
    % Solver only deals with zero homogeneous Dirichlet conditions. Therefore,
    % if nonzero Dirichlet conditions are given, we solve lap(u) = f with u|bc = g
    % as u = v + w, where v|bc = g, and lap(w) = f - lap(v), w|bc = 0:
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
%%

    % Convert rhs to C^{(3/2)} coefficients:
    %F = cheb2ultra( cheb2ultra( F ).' ).';
    
    U = cheb2ultra(U); 
    V = cheb2ultra(V); 
    
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
    U = invDm*U; 
    V = invDn*V; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%  factored Alternating Direction Implicit method %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve TmX + XTn' = F using ADI, which requires O(n^2log(n)log(1/eps))
    % operations:

    % Calculate ADI shifts based on bounds on the eigenvalues of Tn and Tm:
    a = -4/pi^2 * scl_y;
    b = -39*n^-4 * scl_y;
    c = 39*m^-4 * scl_x;
    d = 4/pi^2 * scl_x;
    [p, q] = ADIshifts(a, b, c, d, tol);

    % Run the ADI method:
if use_fiadi==1
    [U, S, V] = fiadi(Tm, -Tn', U, S, V, p, q); 
else
    [U,S, V] = fadi(Tm, -Tn', U, V, p, q); 
end
    %% add compression:
    [U, S, V] = compression(U, diag(S), V, tol); 
    % Convert back to Chebyshev
    U = ultra1mx2cheb(U); 
    V = ultra1mx2cheb(V); 

    
    %% build a chebfun2: 
    U = real(U); %kill any tiny imaginary things
    V = real(V); 
    % better to split S among U and V then to invert and treat as pivots:
    S = diag(S);  
    S = sqrt(S); 
    U = U*diag(S); 
    V = V*diag(S); 
    
    u = chebfun2([]); 
    
    %%trying to directly place the coeffs in the right spot
    % breaks in feval. Not sure why yet.
    %c1techs = chebtech2({' ',U}); 
    %u.cols.funs{1}.onefun = c1techs;
    %c1techs = chebtech2({' ',V});
    %u.rows.funs{1}.onefun = c1techs;
    u.cols = chebfun(U, 'coeffs');      
    u.rows = chebfun(V, 'coeffs');
    u.pivotValues = ones(length(S), 1); 
    u.pivotLocations = nan(length(S), 2); %meaningless here.    
end

end

%%%%%%%%%%%%%%%%%%%%%SUBROUTINES%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
J = ceil( log(16*gam)*log(4/tol)/pi^2 );       % Number of ADI iterations
if ( alp > 1e7 )
    K = (2*log(2)+log(alp)) + (-1+2*log(2)+log(alp))/alp^2/4;
    m1 = 1/alp^2;
    u = (1/2:J-1/2)*K/J;
    dn = sech(u) + .25*m1*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u);
else
    K = ellipke( 1-1/alp^2 );
    [~, ~, dn] = ellipj((1/2:J-1/2)*K/J,1-1/alp^2); % ADI shifts for [-1,-1/t]&[1/t,1]
end
p = T( -alp*dn ); q = T( alp*dn );                  % ADI shifts for [a,b]&[c,d]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% factored ADI: 
function [ZZ, DD, YY] = fadi(A,B,U,V,p,q)

% fadi(A, B, U, V, p, q)
% solves AX - XB = U*V' in low rank form using factored ADI with 
% ADI shift parameters provided by vectors p, q.
% OUTPUT: ZZ*DD*YY' \approx =  X.
% References: 
%
% [1] Benner, Peter, Ren-Cang Li, and Ninoslav Truhar. 
% "On the ADI method for Sylvester equations." J. of Comp. and App. Math.
% 233.4 (2009): 1035-1045. 

[m,r] = size(U);
[n, ~] = size(V);
compute_shifts = 0; 
 
In = speye(n); 
Im = speye(m); 

%check for complex-valued shift parameters
Ns = length(p); 
cp = conj(p); 
cq = conj(q); 

B = B';

%do factored ADI 
Z(:, 1:r) = (A-Im*q(1))\(U);
Y(:,1:r) = (B-cp(1)*In)\V; 
ZZ = Z;  
YY = Y; 
DD =  (q(1)-p(1)).*ones(r,1); 
    for i = 1:Ns-1 
        Z = Z + (A - q(i+1)*Im)\((q(i+1)-p(i))*Z);
        Y = Y + (B-cp(i+1)*In)\((cp(i+1)-cq(i))*Y);
        ZZ = [ZZ Z]; 
        YY = [YY Y];
        DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    end
    DD = diag(DD); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% low rank compression: 
function [A, S, B] = compression(Z, d, Y, tol)
% X = Z*diag(d)*Y';  
% Compresses the factors of X:||X -A*S*B'||_2 \leq tol.

%share out diag if pivots are real
if any(~isreal(d))
    Z= Z*diag(d); 
else
ds = sign(d);
d = sqrt(abs(d));
Z = Z*diag(d); 
Y = (Y*diag(d.*ds));
end

%compression step
[QZ, RZ] = qr( Z, 0); 
[QY, RY] = qr( Y, 0);

% % XX * YY^T = (QX*RX) * (QY*RY)^T = QX*(RX*RY^T)*QY^T
inner_piece = RZ*RY'; 
[A, S, B ] = svd( inner_piece ); 
idx = find(diag(S)>tol, 1, 'last');
A = A(:,1:idx); 
B = B(:,1:idx); 
S = S(1:idx,1:idx); 
A = QZ * A;  
B = QY * B; 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERSION CODES %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = cheb2ultra( X )
% CHEB2ULTRA   Convert vector of Chebyshev coefficients to C^(3/2).
%
%    CHEB2ULTRA(X) applies the conversion to each column of X if X is a
%    matrix.

% First convert the matrix of Chebyshev coefficients to a matrix of
% Legendre coefficients:
m = size( X, 1 );
if ( m <= 10000 ) % Determined experimentally
    S = cheb2leg_mat( m );
    X = S * X;
else
    X = cheb2leg( X );
end

% Now, convert the matrix of Legendre coefficients to a matrix of
% ultraspherical coefficients:
S = leg2ultra_mat( m );
X = S * X;

end

function X = ultra1mx2cheb( X )
% ULTRA1MX2CHEB    Convert vector of (1-x^2)C^(3/2) coefficients to
% Chebyshev.
%
%  ULTRA1MX2CHEB(X) applies the conversion each column of X if X is
%    matrix.

% First, convert the matrix of (1-x^2)C^(3/2)(x) coefficients
% to Legendre coefficients:
m = size( X, 1 );

S = ultra1mx2leg_mat( m );
X = S * X;

% Now, convert the matrix of Legendre coefficient to a matrix of Chebyshev
% coefficients:
if ( m <= 10000 ) % Determined experimentally
    S = leg2cheb_mat( m );
    X = S * X;
else
    X = leg2cheb( X );
end

end

function S = leg2ultra_mat( n )
% LEG2ULTRA_MAT Conversion matrix from Legendre coefficients to C^(3/2).
%
% Given coefficients in the Legendre basis the C^(3/2) coefficients
% can be computed via
%
%     c = rand(10, 1);    % Legendre coefficients
%     S = leg2ultra_mat( length(c) ); % conversion matrix
%     d = S * c;           % C^(3/2) coefficients
%

% Alex Townsend, 5th May 2016

lam = 1/2;
dg = lam./(lam + (2:n-1))';
v  = [1 ; lam./(lam+1) ; dg];
w  = [0 ; 0 ; -dg];
S  = spdiags( [v w], [0 2], n, n );

end

function S = ultra1mx2leg_mat( n )
% ULTRA1MX2LEG_MAT Conversion matrix for (1-x^2)C^(3/2) to Legendre.
%
% Given coefficients in the (1-x^2)C^(3/2) basis the Legendre coefficients
% can be computed via
%
%     c = rand(10, 1);     % (1-x^2)C^(3/2) coefficients
%     S = ultra1mx2leg_mat( length(c) ); % conversion matrix
%     c_leg = S * c;       % Legendre coefficients
%

% Alex Townsend, 5th May 2016

d = ones(n, 1);
S = spdiags(((1:n).*(2:(n+1))./2./(3/2:n+1/2))', 0, n, n);
S = spdiags( [d,-d], [0,-2], n, n ) * S;

end

function L = cheb2leg_mat( N )
% CHEB2LEG_MAT Construct the cheb2leg conversion matrix.

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
% LEG2CHEB_MAT Construct the leg2cheb conversion matrix.

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

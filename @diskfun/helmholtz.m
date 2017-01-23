function u = helmholtz( f, K, bc, m, n )
% HELMHOLTZ   Fast Helmholtz solver with Dirichlet boundary conditions for the disk. 
%   U = HELMHOLTZ(F, K BC, N) solves laplacian(U) + K^2(U) = F on the unit disk, which
%   in polar coordinates (theta, r) is 
%  
%       r^2 u_rr   +   r u_r   +  u_{theta theta} + r^2K^2 u  =   r.^2 * f   
%
%   The equation is discretized on an N x N grid in polar coordinates. The
%   solution is imposed to satisfy the Dirichlet data 
%   U(1,theta) = @(theta) BC(theta). The solution is returned as a 
%   diskfun object. (N must be even.)
%
%   U = HEMLHOLTZ(F,K, BC, M, N) is the same as HELMHOLTZ(F, K,BC, N) but uses a
%   discretization size of M x N. (N must be even.)
%
%   F may be a function handle in polar coordinates, a diskfun, or a matrix 
%   of Fourier-Chebyshev coefficients. 
%  
%
% See also SPHEREFUN/HELMHOLTZ, DISKFUN/POISSON.

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
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if( nargin < 5 )
    n = m;
end
m = 2*m+1; 

% Is n even? 
if ( mod(n, 2) == 1 ) 
    error('DISKFUN:HELMHOLTZ:N', 'The Fourier discretization size must be even.')
end

% Construct operators
D1 = ultraS.diffmat( m, 1 );              % 1st order ultraS diffmat
D2 = ultraS.diffmat( m, 2 );              % 2nd order ultraS diffmat
Mr = ultraS.multmat(m,[0;1],1);           % multiplication of r in ChebU 
Mr2 = ultraS.multmat(m,[.5;0;.5],2);      % multiplication of r^2 in ultra2
Mr2c = ultraS.multmat(m,[.5;0;.5],0);     % multiplication of r^2 in Cheb
S1 = ultraS.convertmat( m, 0, 1 );        % convert chebT coeffs -> ultra2
S12 = ultraS.convertmat( m, 1, 1);        % convert chebU coeffs -> ultra2

% Discretization grid:
x0 = chebpts( m );    
th0 = pi*trigpts( n );

% Forcing term:
if ( isa(f, 'function_handle') )
    f = @(r,th) feval(f, th, r);  % switch convention 
    [rhs_r, rhs_theta] = meshgrid( x0, th0 ); 
    F = rhs_r.^2.*feval( f, rhs_r, rhs_theta );
    % If F is real-valued, we will use symmetry to reduce # computations
    realValued = isreal(F); 
    F = (S1*chebtech2.vals2coeffs( F.' )).';        % Get (C^{(2)},trigvals) basis
    F = trigtech.vals2coeffs( F );
elseif ( isa(f, 'diskfun') )
    realValued = 1; %for now assume diskfuns are real
    F = coeffs2(f, n, m);
    F = S1*Mr2c*F; % r.^2*rhs in C^{2}
    F = F.';  
elseif ( isa( f, 'double' ) ) % assume these are chebyshev coeffs
    %add check if real-valued; for now we will just assume not.
    realValued = 0; 
    f = chebtech2.alias(trigtech.alias(f.', n).',m); % get correct size
    F = S1*Mr2c*f; % r.^2*rhs in C^{2}
    F = F.';
end
 
% If F is real-valued, we will use symmetry to reduce # computations 
% How many Fourier coeffs we need to solve for? n if rv=0, n/2 +1 if rv=1
d = (-n/2)*realValued+n +realValued; 

% Set up LHS
L = Mr2*D2 + S12*Mr*D1 + K^2*S1*Mr2c;
L = L(1:end-2, :); % Eliminate last two rows to make room for BCs.

% Boundary conditions 
bcvals = feval(bc, th0);
bc = trigtech.vals2coeffs(bcvals); 
S1 = S1(1:end-2, :); % Make S1 the right size for including BCs.

% Set up for Sherman-Morrison solve
W = [1 ;zeros((m-1)/2, 1)]; 

CFS = zeros(m,n);
% Main loop
for k = 1 : d
    j = -n/2 + k-1; % Wave number 
    a = mod(j,2)+1; % Even a=1 or odd a=2 wave number
    B= L-j^2*S1 ; 
    
    % We will keep the odd or even contributions only;  
    kB = [zeros(1,(m-1)/2+abs(mod(j,2)-1)); B(a:2:end,a:2:end)];
    kB(1,1)=1;

    % Correct size of vectors (different size system depending on wave
    % number)
    w = W(1:end-(a-1));
    
    % Set up RHS
    b = [bc(k); F(k,a:2:end-2).'];

    % We solve Ax=b by Sherman Morrison formula: 
    % inv(A)b = inv(B + wv^t)b = inv(B)b - ( inv(B)wv^t * inv(B)b)/(1+v^t*inv(B)w)
    invBb = kB\b; 
    invBw = kB\w; 
    
    % We want to compute
    % CFS(a:2:end, k) = invBb-(invBw*v*invBb)/(1+v*invBw); 
    % But this can be made faster by using the fact that v=[0 1 1 1 ... 1]: 
    CFS(a:2:end,k) = invBb-sum(invBb(2:end))*invBw/(1+sum(invBw(2:end)));    
end

% When F is real-valued, this gets the rest of the Fourier coeffs
if ( realValued == 1 )
    CFS(:, d+1:n) = fliplr(conj(CFS(:,2:d-1)));
end

% Solution returned as a diskfun:
u = diskfun.coeffs2diskfun( CFS ); 

end


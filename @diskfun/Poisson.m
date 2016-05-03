
function u = poisson( f, bc, m, n )
% POISSON               FAST POISSON SOLVER FOR THE DISK 
% 
% u = POISSON(F, BC, N) solves  
%  
%       r^2 u_rr   +   r u_r   +  u_{theta theta}   =   r.^2 * f   
%
% on the unit disk written in polar coordinates (r, theta) with Dirichlet
% data  u(1,theta) = @(theta) bc(theta) and a discretization size of N x N.
% The solution is returned as a diskfun object. (N must be even)
%
% u = POISSON(F, BC, M, N) is the same as POISSON(F, BC, N) but uses a
% discretization size of M x N. (N must be even)
%
% EXAMPLE: 
%  bc = @(th) 0*th;              
%  f = @(th, r) -1 + 0*th;            
%  u = diskfun.Poisson( f, bc, 100); 

% DEVELOPER'S NOTE: 
%
% METHOD: Spectral method (in coeff space). We use the Fourier basis in the
% periodic theta-direction and the Ultraspherical spectral method in the
% radial direction.   
%
% LINEAR ALGEBRA: Matrix equations (a generalized Sylvester matrix
% equation). Remarkably, the matrix equation decouples into n linear systems. 
% This are of the form tridiagonal + rank-2 and are solved using the Woodbury
% formula. *Update: using parity properties of Fourier coeffs on the disk,
% this reduces to a tridiagonal + rank-1 (this special case of the Woodbury
% formula is equivalent to the Sherman-Morrison formula)
%
% SOLVE COMPLEXITY:    O( m*n )  N = m*n = total degrees of freedom
%
%   Alex Townsend, July 2015. (Updated Feb 2016,  Heather Wilber)

if ( nargin < 4 ) 
    n = m; 
end

% Double up to use DFS method: 
m = 2*m+1; 
f = @(r,th) feval(f, th, r);  % switch convention 

% Discretization grid:
x0 = chebpts( m );    
th0 = pi*trigpts( n );

[rhs_r, rhs_theta] = meshgrid( x0, th0 ); 
F = rhs_r.^2.*feval( f, rhs_r, rhs_theta );           % Get (chebvals,trigvals) of rhs

%if F is real-valued, we will use symmetry to reduce # computations
rv = isreal(F); 
d = (-n/2)*rv+n +rv; %how many Fourier coeffs we need to solve for? (n if rv=0, n/2 +1 if rv=1)


D1 = ultraS.diffmat( m, 1 );              % 1st order ultraS diffmat
D2 = ultraS.diffmat( m, 2 );              % 2nd order ultraS diffmat
Mr = ultraS.multmat(m,[0;1],1);           % multiplication of r in ChebU 
Mr2 = ultraS.multmat(m,[.5;0;.5],2);      % multiplication of r^2 in ultra2
S1 = ultraS.convertmat( m, 0, 1 );        % convert chebT coeffs -> ultra2
S12 = ultraS.convertmat( m, 1, 1);        % convert chebU coeffs -> ultra2

% set up static elements of LHS
L = Mr2*D2 + S12*Mr*D1;
L = L(1:end-2, :); %eliminate last two rows to make room for BCs.

% Boundary conditions (because of even-odd parity, only one explicit BC is needed) 
bcvals = feval(bc, th0);
bc = trigtech.vals2coeffs(bcvals); 

F = (S1*chebtech2.vals2coeffs( F.' )).';        % Get (C^{(2)},trigvals) basis
F = trigtech.vals2coeffs( F );

S1 = S1(1:end-2, :); %make S1 the right size for including BCs in later computations.

%set up static vectors for Sherman-Morrison solve
W = [1 ;zeros((m-1)/2, 1)]; 



CFS = zeros(n,m);

%main loop
for k = 1 : d
    j = -n/2 + k-1; %wave number 
    a = mod(j,2)+1; %even a=1  or odd a=2 wave number
    B= L-j^2*S1 ; 
    
    %we will keep the odd or even contributions only;  
    kB = [zeros(1,(m-1)/2+abs(mod(j,2)-1)); B(a:2:end,a:2:end)];
    kB(1,1)=1;

    %correct size of vectors (different size system depending
    %on wave number)
    w = W(1:end-(a-1));
    
    
    %set up RHS
    b = [bc(k); F(k,a:2:end-2).'];

    %we solve Ax=b by Sherman Morrison formula: inv(A)b = inv(B + wv^t)b =
    %inv(B)b - ( inv(B)wv^t * inv(B)b)/(1+v^t*inv(B)w)
    
    invBb = kB\b; 
    invBw = kB\w; 
    
    %CFS(k, a:2:end) = invBb-(invBw*v*invBb)/(1+v*invBw); 
    %rewrite using the fact that  v=[0 1 1 1 ... 1]: 
    CFS(k, a:2:end) = invBb-sum(invBb(2:end))*invBw/(1+sum(invBw(2:end)));
    
end


%when F is real-valued, this gets the rest of the Fourier coeffs
if rv==1
CFS(d+1:n,:) = flip(conj(CFS(2:d-1,:)));
end

VALS = trigtech.coeffs2vals( chebtech2.coeffs2vals( CFS.' ).' ); 

% Now restrict down to region of interest:
VALS = VALS(:, (m-1)/2+1:m ); 

% Finally, make a diskfun object out of the values: 
u = diskfun( real( VALS ).' ); 

end


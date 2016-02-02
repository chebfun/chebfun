function u = Poisson( f, bc, m, n )
% POISSON               FAST POISSON SOLVER FOR THE DISK 
% 
% u = POISSON(F, BC, N) solves  
%  
%       r^2 u_rr   +   r u_r   +  u_{theta theta}   =   r.^2 * f   
%
% on the unit disk written in polar coordinates (r, theta) with Dirichlet
% data  u(1,theta) = @(theta) bc(theta) and a discretization size of N x N.
% The solution is returned as a diskfun object. 
%
% u = POISSON(F, BC, M, N) is the same as POISSON(F, BC, N) but uses a
% discretization size of M x N. 
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
% formula. 
%
% SOLVE COMPLEXITY:    O( m*n )  N = m*n = total degrees of freedom
%
%   Alex Townsend, July 2015. 

if ( nargin < 4 ) 
    n = m; 
end

% Double up to use DFS method: 
m = 2*m+1; 
f = @(r,th) feval(f, th, r);  % switch convention 

% Discretization grid:
x0 = chebpts( m );    
th0 = pi*trigpts( n );



% Construct useful spectral matrices: 
DF = -spdiags((-n/2:0)', 0, n/2+1, n/2+1).^2; % 2nd order Fourier diffmat
D1 = ultraS.diffmat( m, 1 );              % 1st order ultraS diffmat
D2 = ultraS.diffmat( m, 2 );              % 2nd order ultraS diffmat
Mr = ultraS.multmat(m,[0;1],1);           % multiplication of r in ChebU 
Mr2 = ultraS.multmat(m,[.5;0;.5],2);      % multiplication of r^2 in ultra2
S1 = ultraS.convertmat( m, 0, 1 );        % convert chebT coeffs -> ultra2
S12 = ultraS.convertmat( m, 1, 1);        % convert chebU coeffs -> ultra2

% Construct discretization for the 'r' part of the disk Laplacian operator: 
L = Mr2*D2 + S12*Mr*D1;

% Are the boundary conditions inhomogeneous or homogenous? 
downvals = feval(bc, th0);
downbc = trigtech.vals2coeffs(downvals); downbc=downbc(1:n/2+1);
% Flip for Dirichlet along r=-1
topbc = trigtech.vals2coeffs([downvals(n/2+1:n) ; downvals(1:n/2)]); 
topbc = topbc(1:n/2+1);

% Forcing term:
% Given a function handle we need to compute the tensor coefficients in the
% (C^{(2)},Fourier) basis. This is a little tricky to get right:
[rhs_r, rhs_theta] = meshgrid( x0, th0 ); 
F = rhs_r.^2.*feval( f, rhs_r, rhs_theta );           % Get (chebvals,trigvals) of rhs
F = (S1*chebtech2.vals2coeffs( F.' )).';        % Get (C^{(2)},trigvals) basis
F = trigtech.vals2coeffs( F ); % fftshift puts them in the  




% Want to solve        X L^T   +  DF X S1^T = F  
%    subject to        X (+/-1's) =  topbc ,   X (1's) = downbc
% Note that the matrix equation decouples because DF is diagonal.
%m = 
% Solve decoupled matrix equation for X, one row at a time: 
CFS = zeros(n/2+1, m); 

%modify BC so even/odd decouples

modified_bc_data1 = (downbc + topbc)/2;
modified_bc_data2 = (downbc - topbc)/2;

for k = 1 : n/2+1
    % Make discretization: 
    A = [                 (-1).^(0:m-1)              ;...
                           ones(1,m)                 ;...
          (L(1:end-2,:) + DF(k,k)*S1(1:end-2,:))    ];
    % Make rhs: 
    
    
    b =  [            topbc(k)      ; ... 
                      downbc(k)     ; ... 
                   F(k,1:end-2).'  ];
    % Solve using Woodbury formula: 
    %CFS(k,:) = WoodburySolver( A, b ); 
    
    % Split into even and odd problems: 
    Aeven = A(1:2:end,1:2:end); 
    Aodd = A(2:2:end, 2:2:end); 
    b =  [     modified_bc_data1(k) ; ... 
               modified_bc_data2(k) ; ... 
                   F(k,1:end-2).'  ];
    % Solve: 
    CFS(k,1:2:end) = WoodburySolver( Aeven, b(1:2:end) );
    CFS(k,2:2:end) = WoodburySolver( Aodd, b(2:2:end) ); 
    
end


CFS=[CFS; flip(conj(CFS(2:end-1,:)))]; %get the rest of the Fourier coeffs

% Convert ot VALUES on the grid: 
VALS = trigtech.coeffs2vals( chebtech2.coeffs2vals( CFS.' ).' ); 

% Now restrict down to region of interest:
VALS = VALS(:, (m-1)/2+1:m ); 

% Finally, make a diskfun object out of the values: 
u = diskfun( real( VALS ).' ); 


% DEBUG 
% % Plot: 
% [rr, tt] = meshgrid( x0, th0 );
% subplot(1,3,1)
% surf(rr, tt, real(VALS),'edgealpha',0,'facecolor','interp'), 
% ylim([-pi pi])
% 
% % Compare to exact solution:
% subplot(1,3,2)
% surf(rr, tt , exact(rr,tt),'edgealpha',0,'facecolor','interp'), 
% ylim([-pi pi])
% 
% subplot(1,3,3)
% surf(rr, tt, abs( VALS - exact( rr, tt) ),'edgealpha',0,'facecolor','interp' )

end

function x = WoodburySolver( A , b )
% WOODBURYSOLVER         FAST SOLVER FOR ALMOST BANDED MATRICES 
% 
% X = WOODBURYSOLVER( A, b) solves Ax=b using Woodbury formula. Here, the
% code is specialized so that A = banded + rank 2, and its the first two
% rows that are dense. 

dense_rows = 2;                            % number of dense rows 
Ik = eye( dense_rows );

% Write A = B + low rank: 
B = A; B(1:dense_rows,:) = 0;        % Make B.  
B(1:dense_rows,1:dense_rows) = Ik; 

% Extract out the dense rows of A: 
V = A(1:dense_rows, :); 

% Low rank part.  A = B + U * V: 
V(1:dense_rows, 1:dense_rows) = V(1:dense_rows, 1:dense_rows) - Ik; 
m = size( A, 1 ); 
U = [ eye( dense_rows ) ; zeros(m-dense_rows, dense_rows) ]; 

% Use the Woodbury formula to solve Ax = b: 
dummy = B \ b; 
x = dummy - B \ ( U * ( ( Ik + V * (B \ U) ) \ (V*dummy) ) );

end 
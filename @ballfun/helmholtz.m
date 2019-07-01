function u = helmholtz(f, K, BC, m, varargin)
%HELMHOLTZ   Helmholtz solver with Dirichlet or Neumann boundary conditions.
%   U = HELMHOLTZ(F, K, BC, M, N, P) is the solution to the Helmholtz
%   equation with right-hand side F, frequency K, and Dirichlet boundary
%   data BC. It uses a discretization size of M x N x P. 
%
%   U = HELMHOLTZ(F, K, BC, M) is the same as HELMHOLTZ(F, K, BC, M, M, M).
%
%   U = HELMHOLTZ(F, K, BC, M, N, P, 'neumann') is the solution to the Helmholtz
%   equation with right-hand side F, frequency K, and Neumann boundary data
%   BC. It uses a discretization size of M x N x P. 
%
%   U = HELMHOLTZ(F, K, BC, M, 'neumann') is the same as 
%   HELMHOLTZ(F, K, BC, M, M, M, 'neumann').
%
%   BC may be a function handle in polar coordinates (lambda, theta), cartesian
%   coordinates (x,y,z), a spherefun, or a matrix of Fourier-Fourier
%   coefficients.
%
% EXAMPLES:
%   f = ballfun(@(r, lam, th) 4 + 4*r.^2.*sin(th).^2,'spherical');
%   exact = @(r, lam, th) r.^2.*sin(th).^2;
%   bc = @(lam, th) exact(1, lam, th);   % Dirichlet condition
%   u = helmholtz(f, 2, bc, 50, 50, 50);
%   exact = ballfun(exact,'spherical');
%   norm(u - exact)
%
%   exact = ballfun(@(x,y,z)(x.^2+y.^2+z.^2).*(x.^2+y.^2));
%   f = laplacian(exact) + 4*exact;
%   bc = spherefun(@(x,y,z)4*(x.^2+y.^2));  % Neumann condition
%   u = helmholtz(f, 2, bc, 40, 'neumann');
%   norm(u - exact)
%
%   Also see POISSON.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPER'S NOTE: 
%
% PROBLEM: We solve:
%
%  r^2sin(th)^2u_rr + 2r(sin(th))^2u_r + (sin(th))^2u_{thth}
%  + cos(th)sin(th)u_th + u_{lamlam} + r^2sin(th)^2K^2 u = (rsin(th)).^2f(r,lam,th)
%
% on the solid sphere (r,lambda, theta) written in spherical coordinates
% with Dirichlet or Neumann condition at r = 1.
%
% METHOD: Spectral method (in coeff space) and Double Fourier sphere method.
% We use the Fourier basis in the theta- and lambda-direction and
% ultraspherical in r.
%
% LINEAR ALGEBRA: Matrix equation solver in (r,theta), using Bartels--Stewart
% algorithm and QZ.
%
% SOLVE COMPLEXITY: O( n^4 )  N = n^3 = total degrees of freedom

% Empty check:
if ( isempty(f) )
    u = f;
    return
end

% Parse user input
isNeumann = any(find(cellfun(@(p) strcmp(p, 'neumann'), varargin)));

if nargin > 4 && isnumeric(varargin{1})
    n = varargin{1};
else
    n = m;
end

if nargin > 5 && isnumeric(varargin{2})
    p = varargin{2};
else
    p = m;
end

if isNeumann
    % Increase m if it's even (Neumann only works for odd m)
    m = m + 1 - mod(m,2);
end

% Get the coeffs
F = coeffs3(f,m,n,p);

% The code was written with variables in the order r, theta, lambda
ord = [1 3 2];
F = permute(F, ord);

% Construct useful spectral matrices (see list above) in r and theta:
DC2 = ultraS.diffmat(m, 2);
S12 = ultraS.convertmat(m, 1, 1);
Mr = ultraS.multmat( m, [0;1], 1);
Mr2 = ultraS.multmat(m, [0.5;0;0.5], 2 );
DC1 = ultraS.diffmat( m, 1);
Msin2 = trigspec.multmat(p, [-.25;0;0.5;0;-0.25] );
I = speye(p);
S02 = ultraS.convertmat(m, 0, 1);
DF2 = trigspec.diffmat(p, 2);
Mcossin = trigspec.multmat(p, [0.25i;0;0;0;-0.25i] );
DF1 = trigspec.diffmat(p, 1);

% Boundary conditions:
[BC1, BC2, bc] = ComputeBoundary(BC, m, n, p, isNeumann);

% Fortunately, the PDE decouples in the lambda variable.
CFS = zeros(m, p, n);

% Solve the Helmholtz equation
if abs(K)>1
    % Divide by K^2
    Lr = Mr2*DC2/K^2 + 2*S12*Mr*DC1/K^2 + Mr2*S02;
else
    Lr = Mr2*DC2 + 2*S12*Mr*DC1 + K^2*Mr2*S02;
end

% Precompute this matrix
Lth = Msin2*DF2 + Mcossin*DF1;

% Normalize bc such that bc(:,2:3) = identity
D = bc(:,2:3);
bc = bc(:,2:3) \ bc;

% Use boundary conditions to remove degrees of freedom 
myS02 = S02;
c1 = myS02(:,2:3);
c2 = Lr(:,2:3);
myS02 = myS02 - myS02(:,2:3)*bc;
Lr = Lr - Lr(:,2:3)*bc;

% Solve the linear system only if f(:,:,k) ~= 0 or BC1(:,k) ~= 0 or BC2(:,k) ~= 0
ListFourierMode = [];
for k = 1:n
    if norm(F(:,:,k),inf) > 1e-16 || norm(BC1(k,:),inf) > 1e-16 || norm(BC2(k,:),inf) > 1e-16
        ListFourierMode = [ListFourierMode k];
    end
end

shift = floor(n/2)+1;

% Loop over the Fourier mode to solve the decoupled equations
for k = ListFourierMode
    
    % Eliminating boundary conditions, changes rhs:
    BC = D \ [BC1(k,:); BC2(k,:)];
    
    % Special case for Poisson equation with Neumann BC
    if( k == floor(n/2)+1 && K == 0 && isNeumann )
        
        % Multiply F by r^2
        ff = Mr2*S02*F(:,:,k);
        
        % Convert the rhs to a Leg x Cheb matrix
        p_tilde = max(2*p-2,1);
        Xchebvals = zeros(p_tilde,m);
        Xchebvals(floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2),:) = ff.';
        Xchebvals = trigtech.coeffs2vals(Xchebvals);
        ff = chebvals2legcoeffs(Xchebvals(1:p,:));
        
        % Convert BC to Leg vector
        Xchebvals = zeros(p_tilde,2);
        Xchebvals(floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2),:) = BC.';
        Xchebvals = trigtech.coeffs2vals(Xchebvals);
        BCLeg = chebvals2legcoeffs(Xchebvals(1:p,:));
               
        % Solution stored in a matrix of Leg x Cheb coeffs
        XLeg = zeros(p,m);
               
        % Loop over the Legendre degrees
        for j = 1:p
            
            % The degrees of freedom from X(:,:,k) have already been removed:
            A = Lr - j*(j-1)*myS02;
            c3 = c2 - j*(j-1)*c1;
            
            % Eliminating boundary conditions, changes rhs:
            ff(j,:) = ff(j,:) - BCLeg(j,:)*c3.';
            
            % Solve for nonzero degree
            if j > 1
                X = A(1:end-2,[1 4:end]) \ ff(j,1:end-2).';
                X = X.';
            % Solve for zero degree
            else
                X = A(1:end-3,4:end) \ ff(j,1:end-3).';
                X = [0 X.'];
            end
            
            % Put the bcs back in:
            col = BCLeg(j,:) - X*bc(:,[1 4:end]).';
            XLeg(j,:) = [X(1), col, X(2:end)];
        end
        
        % Convert back to Fourier x Cheb
        Xchebvals = legcoeffs2chebvals(XLeg);
        Xfourier = trigtech.vals2coeffs([Xchebvals ; flipud(Xchebvals(2:end-1,:))]);
        
        % Fill in the tensor of Fourier x Cheb coeffs
        CFS(:, :, k) = Xfourier(floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2),:).';
        
    else
    
        % Define the operator 
        A = Lth - (k-shift)^2*I;

        % Multiply F by r^2*sin(th)^2
        ff = Mr2*S02*F(:,:,k)*Msin2.';

        if abs(K)>1
            % Divide Lth and f by K^2
            A = A/K^2;
            ff = ff/K^2;
        end

        % Eliminating boundary conditions, changes rhs:
        ff = ff - c1*BC*A.' - c2*BC*Msin2.';

        % Solve resulting Sylvester matrix equation: 
        X = chebop2.bartelsStewart(Lr(1:end-2,[1 4:end]),Msin2,...
            myS02(1:end-2,[1 4:end]),A,ff(1:end-2,:),0,0);

        % Put the bcs back in:
        col = BC - bc(:,[1 4:end])*X;
        CFS(:, :, k) = [X(1,:); col; X(2:end,:)];
    end
end

% Permute back: 
ord = [1 3 2];
CFS = permute( CFS, ord);

% Create ballfun object: 
u = ballfun( CFS, 'coeffs');
end

% Compute boundary rows and Dirichlet and Neumann boundary conditions at
% r = -1
function [BC1, BC2, bc] = ComputeBoundary(BC, m, n, p, isNeumann)
% if BC = function_handle
if isa(BC, 'function_handle')
    % BC = function_handle of x, y and z
    if nargin(BC) == 3
        BC_polar = @(lam,th)BC(cos(lam).*sin(th),sin(lam).*sin(th),cos(th));
        BC = BC_polar;
    end
    % BC is a function handle of lam and theta
    % Grid
    th = pi*trigpts(p);
    lam = pi*trigpts(n);
    % Evaluate function handle at tensor grid:
    [ll, tt] = ndgrid(lam, th);
    BC1 = feval(BC, ll, tt);
    % Test if the function is constant
    if size(BC1) == 1
        BC1 = ones(n,p)*BC1(1);
    end
    % Convert boundary conditions to coeffs:
    BC1 = trigtech.vals2coeffs( trigtech.vals2coeffs( BC1 ).' ).';
% if g is a spherefun
elseif isa(BC, 'spherefun')
    % Get the Fourier x Fourier coefficients of the spherefun
    BC1 = coeffs2(BC, n, p).';
% if g is an array of fourier coefficients lambda x theta
else
    % BC1 is an array of coefficients theta x lambda of size [m,n,p]
    BC1 = trigtech.alias(trigtech.alias(BC.',p).',n);
end

% Dirichlet BC
if ~isNeumann
    % Use the symmetries to find BC2
    BC2 = (-1).^((1:p)-floor(p/2)-1).*BC1;
    
    % Dirichlet boundary rows
    bc1 = ones(1,m);
    bc2 = (-1).^(0:m-1);

% Neumann BC
else
    % BC1 is the derivative of a smooth function on the ball, which contains
    % element of the form r^k exp(i*n*theta) where mod(k,2) = mod(n,2)
    BC2 = (-1).^((1:p)-floor(p/2)).*BC1;
    % Neumann boundary rows
    bc1 = (0:m-1).^2;
    bc2 = (-1).^(1:m).*bc1;
end

% Build boundary conditions matrix
bc = [bc1 ; bc2];
end
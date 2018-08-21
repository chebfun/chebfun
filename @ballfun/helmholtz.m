function u = helmholtz(f, K, g)
%HELMHOLTZ   Helmholtz solver with Dirichlet boundary conditions.
%   HELMHOLTZ(F, K, G) is the solution to the Helmholtz
%   equation with right-hand side F, frequency K, and Dirichlet boundary
%   data given by g(lambda, theta).
%
% Also see POISSON, HELMHOLTZ_NEUMANN.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% PROBLEM: We solve:
%
%  r^2sin(th)^2u_rr + 2r(sin(th))^2u_r + (sin(th))^2u_{thth}
%  + cos(th)sin(th)u_th + u_{lamlam} + r^2sin(th)^2K^2 u = (rsin(th)).^2f(r,lam,th)
%
% on the solid sphere (r,lambda, theta) written in spherical coordinates
% with Dirichlet condition at r = 1.
%
% METHOD: Spectral method (in coeff space) and Double Fourier sphere method.
% We use the Fourier basis in the theta- and lambda-direction and
% ultraspherical in r.
%
% LINEAR ALGEBRA: Matrix equation solver in (r,theta), using Bartels--Stewart
% algorithm and QZ. The matrix decouples in lambda.
%
% SOLVE COMPLEXITY:    O( n^4 )  N = n^3 = total degrees of freedom
%

% DOUBLE FOURIER SPHERE METHOD
%
% Solve:
%
%  r^2sin(th)^2u_rr + 2r(sin(th))^2u_r + (sin(th))^2u_{thth}
%  + cos(th)sin(th)u_th + u_{lamlam} + r^2sin(th)^2K^2 u = r^2sin(th)^2f(r,lam,th)
% on [0,1]x[-pi,pi]x[0,pi],
%
% together with Dirichlet boundary conditions
%
%      u(1,lam,th) = [ g   h ],   (lam,th) in [-pi,pi]x[0,pi].
%
% Using the Double Fourier sphere method, we actually solve:
%
%  r^2sin(th)^2v_rr + 2r(sin(th))^2v_r + (sin(th))^2v_{thth}
%  + cos(th)sin(th)v_th + v_{lamlam} + r^2sin(th)^2K^2 v = (rsin(th)).^2\tilde{f}(r,lam,th)
% on [-1,1]x[-pi,pi]x[-pi,pi],
%
% together with Dirichlet boundary conditions
%
%    u(1,lam,th) = [ g   h ;  flipud(h)  flipud( g )],   (lam,th) in [-pi,pi]x[-pi,pi],
%    u(-1,lam,th) = [ h   g ;  flipud(g)  flipud( h )],  (lam,th) in [-pi,pi]x[-pi,pi].


% DISCRETIZING THE OPERATOR
%
% Laplace operator in spherical coordinates:
% lap =  r^2sin(th)^2u_rr + 2r(sin(th))^2u_r + (sin(th))^2u_{thth}
%                                 + cos(th)sin(th)u_lam + u_{lamlam}
%
% This operator can be discretized as
%
% kron(M_{r^2}D_2^C + 2S_{12}M_rD_1^C, M_{sin(th)^2}, I ) + ...
%             kron(S_{02}, M_{sin(th)^2}D_2^F + M_{cos(th)sin(th)}D_1^F, I )
%             kron(S_{02}, I, D_2^F),
% where
%        M_{r^2} = multmat in ultraS(2) for r^2
%        D_2^C   = 2nd ultraS diffmat
%        S_{12}  = ultraS convertmat from ultraS(1) -> ultraS(2)
%        M_r     = multmat in ultraS(1) for r
%        D_1^C   = 1st ultraS diffmat
%        M_{sin(th)^2} = multmat in Fourier for sin(th)^2
%        I       = identity matrix
%        S_{02}  = ultraS convertmat from chebT -> ultraS(2)
%        D_2^F   = 2nd Fourier diffmat
%        M_{cos(th)sin(th)} = multmat in Fourier for cos(th)sin(th)
%        D_1^F   = 1st Fourier diffmat

[m,n,p] = size(f);

% Get coeffs
F = f.coeffs;

% The code was written with variables in the order theta, r, lambda
ord = [3 1 2];
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
% 2nd Fourier diffmat in lam
DF2lam = trigspec.diffmat(n, 2);

% Boundary conditions:
%
% We know that
%
% v(\pm1,lam,th) = sum_{ijk} c_{ijk} exp(1i*i*lam*pi)*T_j(\pm1)*exp(1i*k*th*pi)
%
%  and
%
% g(lam,th) = sum_{ik} g_{ik} exp(1i*i*lam*pi)*exp(1i*k*th*pi).
%
% Thus, for a fix i and k, we have
%
%         sum_{j} c_{ijk}T_j(1)   =    g_{ik}.

%% Remove this in Julia
% if g = function_handle of lambda, th
if isa(g, 'function_handle')
    % Grid
    th = pi*trigpts(p);
    lam = pi*trigpts(n);
    % Evaluate function handle at tensor grid:
    [ll, tt] = ndgrid(lam, th);
    BC1 = feval(g, ll, tt).';
    % Test if the function is constant
    if size(BC1) == 1
        BC1 = ones(p,n)*BC1(1);
    end
    % Convert boundary conditions to coeffs:
    BC1 = trigtech.vals2coeffs( trigtech.vals2coeffs( BC1 ).' ).';
    
%% if g is an array of fourier coefficients lambda x theta
else
    % BC1 is an array of coefficients theta x lambda
    BC1 = g.';
end

% Use the symmetries to find BC2
BC2 = (-1).^((1:p)-floor(p/2)-1).'.*BC1;

% Fortunately, the PDE decouples in the lambda variable.
CFS = zeros(p, m, n);

% Solve the Helmholtz equation
if abs(K)>1
    % Divide by K^2
    Lr = Mr2*DC2/K^2 + 2*S12*Mr*DC1/K^2 + Mr2*S02;
else
    Lr = Mr2*DC2 + 2*S12*Mr*DC1 + K^2*Mr2*S02;
end

Lth = Msin2*DF2+Mcossin*DF1;
bc1 = (ones(1,m) + (-1).^(0:m-1))/2;
bc2 = (ones(1,m) - (-1).^(0:m-1))/2;
myS02 = S02; 
myLr = Lr; 
c1 = myS02(:,1);
myS02 = myS02 - myS02(:,1)*bc1;
c2 = myS02(:,2);
myS02 = myS02 - myS02(:,2)*bc2;
c3 = myLr(:,1);
myLr = myLr - myLr(:,1)*bc1;
c4 = myLr(:,2);
myLr = myLr - myLr(:,2)*bc2;

%% Solve the linear system only if f(:,:,k) ~= 0 or BC1(:,k) ~= 0 or 
% BC2(:,k) ~= 0
ListFourierMode = [];
for k = 1:n
    if norm(F(:,:,k),inf) > 1e-16 || norm(BC1(:,k),inf) > 1e-16 || norm(BC2(:,k),inf) > 1e-16
        ListFourierMode = [ListFourierMode k];
    end
end

for k = ListFourierMode
    % X(:,:,k) is the solution to 2D PDE:
    %
    %  kron(M_{r^2}D_2^C + 2S_{12}M_rD_1^C, M_{sin(th)^2} ) + ...
    %             kron(S_{02}, M_{sin(th)^2}D_2^F + M_{cos(th)sin(th)}D_1^F + D_2^Flam(k,k)*I )
    
    A = Lth + DF2lam(k,k)*I;
    
    % Multiply F by r^2*sin(th)^2
    ff = Msin2*F(:,:,k)*S02.'*Mr2.';

    if abs(K)>1
        % Divide Lth and f by K^2
        A = A/K^2;
        ff = ff/K^2;
    end
    
    % Eliminating boundary conditions, changes rhs:
    ff = ff - (A*(BC1(:,k)+BC2(:,k))/2)*c1';
    ff = ff - (A*(BC1(:,k)-BC2(:,k))/2)*c2';
    ff = ff - (Msin2*(BC1(:,k)+BC2(:,k))/2)*c3';
    ff = ff - (Msin2*(BC1(:,k)-BC2(:,k))/2)*c4';
    
    % Solve resulting Sylvester matrix equation: 
    X = chebop2.bartelsStewart(Msin2,myLr(1:end-2,3:end),A,...
        myS02(1:end-2,3:end),ff(:,1:end-2),0,0);
    
    % Put the bcs back in:
    col1 =  (BC1(:,k)+BC2(:,k))/2 - X * bc1(3:end).';
    col2 =  (BC1(:,k)-BC2(:,k))/2 - X * bc2(3:end).';
    
    CFS(:, :, k) = [col1 col2 X];
end

% Permute back: 
ord = [2 3 1];
CFS = permute( CFS, ord);

% Create ballfun object: 
u = ballfun( CFS ); 
end

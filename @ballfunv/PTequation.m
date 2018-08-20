function u = PTequation(f) 
% PTEQUATION Solve the equations to compute the poloidal-toroidal 
% decompositon of a BALLFUNV
%   PTEQUATION(f, varargin) is the solution to the equation 
%   sin^2(theta)Delta_S u = f with the condition that the 0th Fourier mode  
%   in lambda and theta is 0. 
 
% Solve the equation with optimal complexity O(m2n^3), where n is the 
% level of discretization in one variable and m = 2 is the Fourier number 
% of sin(th)^2 and sin(th)*cos(th). m corresponds to the bandwidth of  
% A = sin^2(theta)Delta_S  
 
[m,n,p] = size(f); 

% Get coeffs 
F = f.coeffs; 
 
% The code was written with variables in the order theta, r, lambda 
ord = [3 1 2]; 
F = permute(F, ord);
 
% Construct useful spectral matrices (see list above): 
Msin2 = trigspec.multmat(p, [-.25;0;.5;0;-.25] ); 
DF1 = 1i*spdiags((-floor(p/2):floor(p/2))', 0, p, p); 
DF2 = trigspec.diffmat(p, 2); 
Mcossin = trigspec.multmat(p, [.25i;0;0;0;-.25i] ); 
I = speye(p,p); 
DF2lam = trigspec.diffmat(n, 2); 
 
% Fortunately, the PDE decouples in the lambda variable 
Lth = Msin2*DF2+Mcossin*DF1; 
 
% Coefficients of the solution 
CFS = zeros(p, m, n); 
 
% Solve the linear system only if f(:,:,k) ~= 0  
ListFourierMode = []; 
for k = 1:n 
    if norm(F(:,:,k),inf) > 1e-16 
        ListFourierMode = [ListFourierMode k]; 
    end 
end 
 
% Decouple in lambda 
for k = ListFourierMode 
     
    % Compute the right hand side term 
    ff = F(:,:,k); 
    
    % Linear system 
    Lthlam = Lth+DF2lam(k,k)*I; 
     
    % Add the row condition to the 0th Fourier mode in lambda 
    if k == floor(n/2)+1 
        Lthlam(floor(p/2)+1,:) = zeros(1,p); 
        Lthlam(floor(p/2)+1,floor(p/2)+1) = 1; 
        % Add the row conditions to rhs 
        ff(floor(p/2)+1,:) = zeros(1,m); 
    end 
     
    % Solve the linear system only if f(:,i,k) ~= 0  
    ListChebyshev = []; 
    for i = 1:m 
        if norm(ff(:,i),inf) > 1e-16 
            ListChebyshev = [ListChebyshev i]; 
        end 
    end 
     
    % Solve the system 
    CFS(:,ListChebyshev,k) = Lthlam\ff(:,ListChebyshev); 
end 
 
% Permute back 
ord = [2 3 1]; 
CFS = permute(CFS, ord); 
 
% Return the solution 
u = ballfun(CFS); 
end 

function varargout = PTdecomposition(v)
% PTDECOMPOSITION Poloidal-toroidal decomposition of a BALLFUNV
%   [P,T] = PTDECOMPOSITION(V) returns the poloidal and toroidal BALLFUN of 
%   the BALLFUNV V, where V = curl(curl(rP)) + curl(rT).
%   r denotes the vector of length sqrt(x^2+y^2+z^2) in the radial
%   direction.
% 
%   Also see PT2BALLFUNV.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check if one component is empty
if isempty(v)
    isVempty = 1;
else
    isVempty = cellfun(@isempty, v.comp, 'UniformOutput', false);
    isVempty = any(cell2mat(isVempty));
end
if isVempty
    error('BALLFUNV:PTdecomposition:input', ...
          'ballfunv must not have an empty component');    
end

% This slows down the code a lot
if norm(div(v)) > 1e-8
    warning('CHEBFUN:BALLFUN:PTdecomposition: the vector field is not divergence-free');
end

% Get the discretization : take the maximum over the components of v
S = max(size(v),[],1);
m = S(1)+1; n = S(2)+2; p = S(3)+6;

% Spectral matrices
Mr = ultraS.multmat(m, [0;1], 0);
MsinL = trigspec.multmat(n, [0.5i;0;-0.5i]); 
McosL = trigspec.multmat(n, [0.5;0;0.5]);
DF1L = 1i*spdiags((-floor(n/2):floor(n/2))', 0, n,n);

MsinT = trigspec.multmat(p, [0.5i;0;-0.5i]);
Msin2T = trigspec.multmat(p, [-0.25;0;0.5;0;-0.25]);
McossinT = trigspec.multmat(p, [0.25i;0;0;0;-0.25i]);
Msin2cosT = trigspec.multmat(p, [-.125;0;.125;0;.125;0;-.125]);
Msin3T = trigspec.multmat(p, [-.125i;0;.375i;0;-.375i;0;.125i]);
DF1T = 1i*spdiags((-floor(p/2):floor(p/2))', 0, p, p);

% Extract coeffs
[Vx,Vy,Vz] = v.comp{:};
Vx = coeffs3(Vx,m,n,p);
Vy = coeffs3(Vy,m,n,p);
Vz = coeffs3(Vz,m,n,p);

% Permute Vx, Vy and Vz
Vx = permute(Vx,[2,3,1]);
Vy = permute(Vy,[2,3,1]);
Vz = permute(Vz,[2,3,1]);

% Rhs of the poloidal part
RhsP = zeros(size(Vx));
% Rhs of the toroidal part
RhsT = zeros(size(Vx));

% Loop over r
for k = 1:m
   RhsP(:,:,k) = -McosL*Vx(:,:,k)*Msin3T.' - MsinL*Vy(:,:,k)*Msin3T.' - Vz(:,:,k)*Msin2cosT.';
   RhsT(:,:,k) = MsinL*Vx(:,:,k)*MsinT.'*DF1T.'*MsinT.' + DF1L*McosL*Vx(:,:,k)*McossinT.'...
               - McosL*Vy(:,:,k)*MsinT.'*DF1T.'*MsinT.' + DF1L*MsinL*Vy(:,:,k)*McossinT.'...
               - DF1L*Vz(:,:,k)*Msin2T.';
end

% Permute back
RhsP = permute(RhsP,[3,1,2]);
RhsT = permute(RhsT,[3,1,2]);

% Multiply RhsP by r
for k = 1:p
   RhsP(:,:,k) = Mr*RhsP(:,:,k);
end

% Poloidal and toroidal scalars
P = PTequation(ballfun(RhsP,'coeffs'));
T = PTequation(ballfun(RhsT,'coeffs'));

% Prepare output:
if ( nargout <= 1 ) 
    varargout = { [P,T] };
else 
    varargout = { P, T };
end
end

function u = PTequation(f) 
% PTEQUATION Solve the equations to compute the poloidal-toroidal 
% decompositon of a BALLFUNV
%   PTEQUATION(f, varargin) is the solution to the equation 
%   sin^2(theta)Delta_S u = f with the condition that the 0th Fourier mode  
%   in lambda and theta is 0. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
u = ballfun(CFS,'coeffs'); 
end 


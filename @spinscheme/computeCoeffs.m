function coeffs = computeCoeffs(K, dt, L, LR, S)
%COMPUTECOEFFS   Compute coefficients of a SPINSCHEME.
%   COEFFS = COMPUTECOEFFS(K, DT, L, LR, S) computes the coefficients needed by  
%   the SPINSCHEME K from the timestep DT, the linear part L, the linear part 
%   for complex means LR, and the SPINOP S.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
s = K.internalStages;
q = K.steps;
dim = S.dimension;
nVars = S.numVars;
schemeName = K.scheme;
A = cell(s);
B = cell(s, 1);
C = zeros(s, 1);
U = cell(s, q-1);
V = cell(q-1, 1);
phit = cell(s);
N = size(L, 1)/nVars;

% Compute the first three phi-functions (needed by all schemes):
phi1 = spinscheme.phiEval(1, LR, N, dim, nVars);
phi2 = spinscheme.phiEval(2, LR, N, dim, nVars);
phi3 = spinscheme.phiEval(3, LR, N, dim, nVars);

% Take real part for diffusive problems (real eigenvalues):
if ( isreal(L) == 1 )
    phi1 = real(phi1);
    phi2 = real(phi2);
    phi3 = real(phi3);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETD MULTISTEP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETD RUNGE-KUTTA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmpi(schemeName, 'etdrk4') == 1 )
        
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = phit{1,2};
    A{4,3} = 2*phit{1,2};
    
    % Compute B:
    B{2} = 2*phi2 - 4*phi3;
    B{3} = 2*phi2 - 4*phi3;
    B{4} = -phi2 + 4*phi3;

elseif ( strcmpi(schemeName, 'exprk5s8') == 1 ) 
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1/4;
    C(5) = 1/2;
    C(6) = 1/5;
    C(7) = 2/3;
    C(8) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = phit{1,2};
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{1,5} = phit{1,2};
    phit{1,6} = spinscheme.phitEval(1, C(6), LR, N, dim, nVars);
    phit{1,7} = spinscheme.phitEval(1, C(7), LR, N, dim, nVars);
    phit{1,8} = phi1;
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{2,4} = spinscheme.phitEval(2, C(4), LR, N, dim, nVars);
    phit{2,6} = spinscheme.phitEval(2, C(6), LR, N, dim, nVars);
    phit{2,7} = spinscheme.phitEval(2, C(7), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);
    phit{3,6} = spinscheme.phitEval(3, C(6), LR, N, dim, nVars);
    phit{3,7} = spinscheme.phitEval(3, C(7), LR, N, dim, nVars);
    phit{4,6} = spinscheme.phitEval(4, C(6), LR, N, dim, nVars);
    phit{4,7} = spinscheme.phitEval(4, C(7), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*phit{2,2};
    A{4,3} = 2*phit{2,4};
    A{5,3} = -2*phit{2,2} + 16*phit{3,2};
    A{5,4} = 8*phit{2,2} - 32*phit{3,2};
    A{6,4} = 8*phit{2,6} - 32*phit{3,6};
    A{7,4} = -(125/162)*A{6,4};
    A{6,5} = -2*phit{2,6} + 16*phit{3,6};
    A{7,5} = (125/1944)*A{6,4} - (4/3)*phit{2,7} + (40/3)*phit{3,7};
    Phi = (5/32)*A{6,4} - (25/28)*phit{2,6} + (81/175)*phit{2,7} - ...
        (162/25)*phit{3,7} + (150/7)*phit{4,6} + (972/35)*phit{4,7} + 6*phi4;
    A{8,5} = -(16/3)*phi2 + (208/3)*phi3 - 40*Phi;
    A{7,6} = (3125/3888)*A{6,4} + (25/3)*phit{2,7} - (100/3)*phit{3,7};
    A{8,6} = (250/21)*phi2 - (250/3)*phi3 + (250/7)*Phi;
    A{8,7} = (27/14)*phi2 - 27*phi3 + (135/7)*Phi;
    
    % Compute B:
    B{6} = (125/14)*phi2 - (625/14)*phi3 + (1125/14)*phi4;
    B{7} = -(27/14)*phi2 + (162/7)*phi3 - (405/7)*phi4;
    B{8} = (1/2)*phi2 - (13/2)*phi3 + (45/2)*phi4;
    
elseif ( strcmpi(schemeName, 'friedli') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*phit{2,2};
    A{4,2} = -26/25*phi1 + 2/25*phi2;
    A{4,3} = 26/25*phi1 + 48/25*phi2;
    
    % Compute B:
    B{1} = phi1 - 3*phi2 + 4*phi3;
    B{3} = 4*phi2 - 8*phi3;
    B{4} = -phi2 + 4*phi3;
    
elseif ( strcmpi(schemeName, 'hochbruck-ostermann') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    C(5) = 1/2;
    
    % Compute the phi functions:
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{1,5} = spinscheme.phitEval(1, C(5), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);
     
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4*phit{2,2};
    A{4,2} = phi2;
    A{5,2} = 1/4*phi2 - phi3 + 2*phit{2,2} - 4*phit{3,2};
    A{4,3} = A{4,2};
    A{5,3} = A{5,2};
    A{5,4} = phit{2,2} - A{5,2};
    
    % Compute B:
    B{1} = phi1 - 3*phi2 + 4*phi3;
    B{4} = -phi2 + 4*phi3;
    B{5} = 4*phi2 - 8*phi3;
    
elseif ( strcmpi(schemeName, 'krogstad') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4*phit{2,2};
    A{4,3} = 2*phi2;
    
    % Compute B:
    B{2} = 2*phi2 - 4*phi3;
    B{3} = 2*phi2 - 4*phi3;
    B{4} = -phi2 + 4*phi3;
    
elseif ( strcmpi(schemeName, 'minchev') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4/25*phit{1,2} + 24/25*phit{2,2};
    A{4,2} = 21/5*phi2 - 108/5*phi3;
    A{4,3} = 1/20*phi1 - 33/10*phi2 + 123/5*phi3;
    
    % Compute B:
    B{1} = 31/30*phi1 - 17/5*phi2 + 6*phi3 - 4*phi4;
    B{2} = -1/10*phi1 + 1/5*phi2 - 4*phi3 + 12*phi4;
    B{3} = 1/30*phi1 + 23/5*phi2 - 8*phi3 - 4*phi4;
    B{4} = 1/30*phi1 - 7/5*phi2 + 6*phi3 - 4*phi4;
    
elseif ( strcmpi(schemeName, 'strehmel-weiner') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*phit{2,2};
    A{4,2} = -2*phi2;
    A{4,3} = 4*phi2;
    
    % Compute B:
    B{1} = phi1 - 3*phi2 + 4*phi3;
    B{3} = 4*phi2 - 8*phi3;
    B{4} = -phi2 + 4*phi3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPI RUNGE-KUTTA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAWSON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'lawson4') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi0 = spinscheme.phiEval(0, LR, N, dim, nVars);
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi0 = real(phi0);
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{1} = 1/6*phi0;
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERALIZED LAWSON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'genlawson41') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi0 = spinscheme.phiEval(0, LR, N, dim, nVars);
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);

    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi0 = real(phi0);
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
       
elseif ( strcmpi(schemeName, 'genlawson42') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -phit{2,2};
    U{3,1} = -phit{2,2} + 1/4;
    U{4,1} = -phi2 + 1/2*phit02;

    % Compute V:
    V{1} = -phi2 + 1/3*phit02 + 1/6;
    
elseif ( strcmpi(schemeName, 'genlawson43') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);

    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -2*phit{2,2} - 2*phit{3,2};
    U{3,1} = -2*phit{2,2} - 2*phit{3,2} + 5/8;
    U{4,1} = -2*phi2 - 2*phi3 + 5/4*phit02;
    U{2,2} = 1/2*phit{2,2} + phit{3,2};
    U{3,2} = 1/2*phit{2,2} + phit{3,2} - 3/16;
    U{4,2} = 1/2*phi2 + phi3 - 3/8*phit02;
    
    % Compute V:
    V{1} = -2*phi2 - 2*phi3 + 5/6*phit02 + 1/2;
    V{2} = 1/2*phi2 + phi3 - 1/4*phit02 - 1/6;
      
elseif ( strcmpi(schemeName, 'genlawson44') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);
    phit{4,2} = spinscheme.phitEval(4, C(2), LR, N, dim, nVars);

    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -3*phit{2,2} - 5*phit{3,2} - 3*phit{4,2};
    U{2,2} = 3/2*phit{2,2} + 4*phit{3,2} + 3*phit{4,2};
    U{2,3} = -1/3*phit{2,2} - 1*phit{3,2} - 1*phit{4,2};
    U{3,1} = U{2,1} + 35/32;
    U{3,2} = U{2,2} - 21/32;
    U{3,3} = U{2,3} + 5/32;
    U{4,1} = -3*phi2 - 5*phi3 - 3*phi4 + 35/16*phit02;
    U{4,2} = 3/2*phi2 + 4*phi3 + 3*phi4 - 21/16*phit02;
    U{4,3} = -1/3*phi2 - phi3 - phi4 + 5/16*phit02;
      
    % Compute V:
    V{1} = -3*phi2 - 5*phi3 - 3*phi4 + 35/24*phit02 + 1;
    V{2} = 3/2*phi2 + 4*phi3 + 3*phi4 - 7/8*phit02 - 2/3;
    V{3} = -1/3*phi2 - 1*phi3 - 1*phi4 + 5/24*phit02 + 1/6;
      
        
elseif ( strcmpi(schemeName, 'genlawson45') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phi5 = spinscheme.phiEval(5, LR, N, dim, nVars);
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);
    phit{4,2} = spinscheme.phitEval(4, C(2), LR, N, dim, nVars);
    phit{5,2} = spinscheme.phitEval(4, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phi5 = real(phi5);
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -4*phit{2,2} - 13/12*8*phit{3,2} - 9*phit{4,2} - 4*phit{5,2};
    U{2,2} = 3*phit{2,2} + 19/2*phit{3,2} + 12*phit{4,2} + 6*phit{5,2};
    U{2,3} = -4/3*phit{2,2} - 7/12*16*phit{3,2} - 7*phit{4,2} - 4*phit{5,2};
    U{2,4} = 1/4*phit{2,2} + 88/96*phit{3,2} + 3/2*phit{4,2} + phit{5,2};
    
    U{3,1} = U{2,1} + 105/64;
    U{3,2} = U{2,2} - 189/128;
    U{3,3} = U{2,3} + 45/64;
    U{3,4} = U{2,4} - 35/256;
    
    U{4,1} = -4*phi2 - 26/3*phi3 - 9*phi4 - 4*phi5 + 105/32*phit02;
    U{4,2} = 3*phi2 + 19/2*phi3 + 12*phi4 + 6*phi5 - 189/64*phit02;
    U{4,3} = -4/3*phi2 - 14/3*phi3 - 7*phi4 - 4*phi5 + 45/32*phit02;
    U{4,4} = 1/4*phi2 + 11/12*phi3 + 3/2*phi4 + phi5 - 35/128*phit02;
    
%     u = { one, [], [], [], []; ...
%       ez2, -phi_22 - 13/12*phi_32 - 9/16*phi_42 - 1/8*phi_52, ...
%       3/4*phi_22 + 19/16*phi_32 + 3/4*phi_42 + 3/16*phi_52, ...
%      -1/3*phi_22 - 7/12*phi_32 - 7/16*phi_42 - 1/8*phi_52, ...
%       1/16*phi_22 + 11/96*phi_32 + 3/32*phi_42 + 1/32*phi_52; ...
%       
%       ez2, -phi_22 - 13/12*phi_32 - 9/16*phi_42 - 1/8*phi_52 + 105/64*one, ...
%       3/4*phi_22 + 19/16*phi_32 + 3/4*phi_42 + 3/16*phi_52 - 189/128*one, ...
%      -1/3*phi_22 - 7/12*phi_32 - 7/16*phi_42 - 1/8*phi_52 + 45/64*one, ...
%       1/16*phi_22 + 11/96*phi_32 + 3/32*phi_42 + 1/32*phi_52 - 35/256*one; ...
%       
%       ez, -4*phi_2 - 26/3*phi_3 - 9*phi_4 - 4*phi_5 + 105/32*ez2, ...
%       3*phi_2 + 19/2*phi_3 + 12*phi_4 + 6*phi_5 - 189/64*ez2, ...
%      -4/3*phi_2 - 14/3*phi_3 - 7*phi_4 - 4*phi_5 + 45/32*ez2, ...
%       1/4*phi_2 + 11/12*phi_3 + 3/2*phi_4 + phi_5 - 35/128*ez2 };
  
    % Compute V:
    V{1} = -4*phi2 - 26/3*phi3 - 9*phi4 - 4*phi5 + 35/16*phit02 + 5/3;
    V{2} = 3*phi2 + 19/2*phi3 + 12*phi4 + 6*phi5 - 63/32*phit02 - 5/3;
    V{3} = -4/3*phi2 - 14/3*phi3 - 7*phi4 - 4*phi5 + 15/16*phit02 + 5/6;
    V{4} = 1/4*phi2 + 11/12*phi3 + 3/2*phi4 + phi5 - 35/192*phit02 - 1/6;
    
%     v = { ez,  -4*phi_2 - 26/3*phi_3 - 9*phi_4 - 4*phi_5 + 5/3*one + 35/16*ez2, ...
%       3*phi_2 + 19/2*phi_3 + 12*phi_4 + 6*phi_5 - 5/3*one - 63/32*ez2, ... 	
%      -4/3*phi_2 - 14/3*phi_3 - 7*phi_4 - 4*phi_5 + 5/6*one + 15/16*ez2, ...
%       1/4*phi_2 + 11/12*phi_3 + 3/2*phi_4 + phi_5 - 1/6*one - 35/192*ez2; ...
%       [], [], [], [], []; ...
%       [], one, [], [], []; ...
%       [], [], one, [], []; ...
%       [], [], [], one, [] };
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFIED GENERALIZED LAWSON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'modgenlawson43') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phit02 = spinscheme.phitEval(0, C(2), LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{1,4} = spinscheme.phitEval(1, C(4), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);

    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit02 = real(phit02);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = phit02;
    
    % Compute B:
    B{2} = 1/3*phit02;
    B{3} = B{2};
    B{4} = 1/3*phi2 + phi3 + phi4 - 5/24*phit02;
    
    % Compute U:
    U{2,1} = -2*phit{2,2} - 2*phit{3,2};
    U{3,1} = -2*phit{2,2} - 2*phit{3,2} + 5/8;
    U{4,1} = -2*phi2 - 2*phi3 + 5/4*phit02;
    U{2,2} = 1/2*phit{2,2} + phit{3,2};
    U{3,2} = 1/2*phit{2,2} + phit{3,2} - 3/16;
    U{4,2} = 1/2*phi2 + phi3 - 3/8*phit02;
    
    % Compute V:
    V{1} = -phi2 + phi3 + 3*phi4 + 5/24*phit02;
    V{2} = 1/6*phi2 - phi4 - 1/24*phit02;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREDICTOR-CORRECTOR:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif ( strcmpi(schemeName, 'pecec433') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi- and phit-functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/3*phi2 + phi3 + phi4;
    
    % Compute B:
    B{3} = 1/3*phi2 + phi3 + phi4;
    
    % Compute U:
    U{2,1} = -2*phi2 - 2*phi3;
    U{3,1} = -phi2 + phi3 + 3*phi4;
    U{2,2} = 1/2*phi2 + phi3;
    U{3,2} = 1/6*phi2 - phi4;
    
    % Compute V:
    V{1} = -phi2 + phi3 + 3*phi4;
    V{2} = 1/6*phi2 - phi4;
    
elseif ( strcmpi(schemeName, 'pecec635') == 1 )
        
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi- and phit-functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phi5 = spinscheme.phiEval(5, LR, N, dim, nVars);
    phi6 = spinscheme.phiEval(6, LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/5*phi2 + 5/6*phi3 + 7/4*phi4 + 2*phi5 + phi6;

    % Compute B:
    B{1} =  phi1 + 13/12*phi2 - 5/4*phi3 - 25/4*phi4 - 9*phi5 - 5*phi6;
    B{3} = 1/5*phi2 + 5/6*phi3 + 7/4*phi4 + 2*phi5 + phi6;
   
    % Compute U:
    U{2,1} = -4*phi2 - 26/3*phi3 - 9*phi4 - 4*phi5;
    U{2,2} = 3*phi2 + 19/2*phi3 + 12*phi4 + 6*phi5;
    U{2,3} = -4/3*phi2 - 14/3*phi3 - 7*phi4 - 4*phi5;
    U{2,4} = 1/4*phi2 + 11/12*phi3 + 3/2*phi4 + phi5;
    U{3,1} = -2*phi2 - 1/3*phi3 + 17/2*phi4 + 16*phi5 + 10*phi6;
    U{3,2} = phi2 + 7/6*phi3 - 11/2*phi4 - 14*phi5 - 10*phi6;
    U{3,3} = -1/3*phi2 - 1/2*phi3 + 7/4*phi4 + 6*phi5 + 5*phi6;
    U{3,4} = 1/20*phi2 + 1/12*phi3 - 1/4*phi4 - phi5 - phi6;
    
    % Compute V:
    V{1} = -2*phi2 - 1/3*phi3 + 17/2*phi4 + 16*phi5 + 10*phi6;
    V{2} = phi2 + 7/6*phi3 - 11/2*phi4 - 14*phi5 - 10*phi6;
    V{3} = -1/3*phi2 - 1/2*phi3 + 7/4*phi4 + 6*phi5 + 5*phi6;
    V{4} = 1/20*phi2 + 1/12*phi3 - 1/4*phi4 - phi5 - phi6;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MISC:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif ( strcmpi(schemeName, 'eglm433') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1;
    
    % Compute the phi functions:
    phi4 = spinscheme.phiEval(4, LR, N, dim, nVars);
    phi5 = spinscheme.phiEval(5, LR, N, dim, nVars);
    phit{1,2} = spinscheme.phitEval(1, C(2), LR, N, dim, nVars);
    phit{1,3} = spinscheme.phitEval(1, C(3), LR, N, dim, nVars);
    phit{2,2} = spinscheme.phitEval(2, C(2), LR, N, dim, nVars);
    phit{3,2} = spinscheme.phitEval(3, C(2), LR, N, dim, nVars);
    
    % Take real part of diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi4 = real(phi4);
        phi5 = real(phi5);
        phit = cellfun(@(f) real(f), phit, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 16/15*phi2 + 16/5*phi3 + 16/5*phi4;
    
    % Compute B:
    B{2} = 32/15*(phi2 + phi3) - 64/5*phi4 - 128/5*phi5;
    B{3} = -1/3*phi2 + 1/3*phi3 + 5*phi4 + 8*phi5;
    
    % Compute U:
    U{2,1} = -2*(phit{2,2} + phit{3,2});
    U{3,1} = -2/3*phi2 + 2*phi3 + 4*phi4;
    U{2,2} = 1/2*phit{2,2} + phit{3,2};
    U{3,2} = 1/10*phi2 - 1/5*phi3 - 6/5*phi4;
    
    % Compute V:
    V{1} = -1/3*phi2 + 5/3*phi3 - phi4 - 8*phi5;
    V{2} = 1/30*phi2 - 2/15*phi3 - 1/5*phi4 + 8/5*phi5;
    
end

% PUT everything in COEFFS:
coeffs.A = A;
coeffs.B = B;
coeffs.C = C;
coeffs.U = U;
coeffs.V = V;

% Compute the missing oefficients using the summation properties of the coeffs:
coeffs = computeMissingCoeffs(K, L, coeffs, dt, phi1, phit);

end

function coeffs = computeMissingCoeffs(K, L, coeffs, dt, phi1, phit)
%COMPUTEMISSINGCOEFFS   Compute the missing oefficients of a SPINSCHEME using 
%the summation properties of the coefficients.
%   COEFFS = COMPUTEMISSINGCOEFFS(K, L, COEFFS, DT, PHI1, PHIT) uses the row 
%   summation properties to compute the COEFFS.A{I,1}, COEFFS.B{1} and COEFFS.E
%   coefficients of the SPINSCHEME K, using the linear part L of the opeartor, 
%   the PHI1 and PHIT functions, and the timestep DT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the coefficients:
A = coeffs.A;
B = coeffs.B; 
C = coeffs.C;
U = coeffs.U;
V = coeffs.V;

% Number of internal stages S and number of steps used Q:
s = K.internalStages;
q = K.steps;

% Precompute the coefficients Ai1 using the row summing property.
for i = 2:s
    A{i,1} = phit{1,i};
    for j = 2:i-1
        if ( ~isempty(A{i,j}) )
            A{i,1} = A{i,1} - A{i,j};
        end
    end
    for j = 1:q-1
        if ( ~isempty(U{i,j}) )
            A{i,1} = A{i,1} - U{i,j};
        end
    end
end

% Precompute the coefficient B1 using the row summing property.
B{1} = phi1;
for i = 2:s
    if ( ~isempty(B{i}) )
        B{1} = B{1} - B{i};
    end
end
for i = 1:q-1
    if ( ~isempty(V{i}) )
        B{1} = B{1} - V{i};
    end
end

% Precompute the E quantities.
E = cell(s+1, 1);
for i = 1:s
   E{i} = exp(C(i)*dt*L);
end
E{s+1} = exp(dt*L);

% Multiply by timestep:
A = cellfun(@(A) dt*A, A, 'UniformOutput', 0);
B = cellfun(@(B) dt*B, B, 'UniformOutput', 0);
U = cellfun(@(U) dt*U, U, 'UniformOutput', 0);
V = cellfun(@(V) dt*V, V, 'UniformOutput', 0);

% Put everything in COEFFS:
coeffs.A = A;
coeffs.B = B;
coeffs.E = E;
coeffs.U = U;
coeffs.V = V;

end
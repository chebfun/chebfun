function schemeCoeffs = computeCoeffs(K, dt, L, M, S)
%COMPUTECOEFFS   Compute coefficients of an EXPINTEG.
%   SCHEMECOEFFS = COMPUTECOEFFS(K, DT, L, M, S) computes the coefficients
%   needed by the EXPINTEG K from the time-step DT, the linear part L, the
%   number of points for complex means M, and the SPINOPERATOR S.
%
% See also EXPINTEG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note (1): At the end of the method the coefficients will be stored in 
% SCHEMECOEFFS, a STRUCT object with the following fields:
%
% SCHEMECOEFFS.A -> sxs CELL-ARRAY that stores the coefficients A, functions of 
%                   the phi-functions
% SCHEMECOEFFS.B -> sx1 CELL-ARRAY that stores the coefficients B, functions of
%                   the phi-functions
% SCHEMECOEFFS.C -> sx1 CELL-ARRAY that stores the numbers C
% SCHEMECOEFFS.E -> (s+1)x1 CELL-ARRAY that stores the s matrix exponentials
%                    e^{C(i)*dt*L}, i=1,..,s, and e^{*dt*L}
% SCHEMECOEFFS.U -> sx(q-1) CELL-ARRAY that stores the coefficients U, functions 
%                   of the phi-functions
% SCHEMECOEFFS.V -> (q-1)x1 CELL-ARRAY that stores the coefficients V, functions 
%                   of the phi-functions
%
% where s is the number of number of internal stages and q is number of steps 
% (>1 for multistep methods). For informations about what these coefficients
% are, see [1].
%
% [1] H. Montanelli and N. Bootland, Solving periodic semilinear stiff PDEs in 
% 1D/2D/3D with exponential integrators, submitted (2017). 

% Note (2): EXPINTEG schemes are used for solving PDEs with diagonal linear 
% operators, e.g., periodic PDEs in 1D/2D/3D with SPIN/SPIN2/SPIN3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE-PROCESSING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set-up:
s = K.stages;          % number of stages
q = K.steps;           % number of steps (>1 for multistep methods)
dim = getDimension(S); % spatial dimension (1, 2 or 3)
nVars = S.numVars;     % number of unknown functions
schemeName = K.scheme; % scheme
N = size(L, 1)/nVars;  % grid points

% Initialize the coefficients of the scheme:
A = cell(s);
B = cell(s, 1);
C = zeros(s, 1);
U = cell(s, q-1);
V = cell(q-1, 1);

% Create a contour around each eigenvalue of the linear part L:
LR = computeLR(S, dt, L, M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETD ADAMS-BASHFORTH:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmpi(schemeName, 'abnorsett4') == 1 )
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,1} = expinteg.psiEval(1, C(1), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute V:
    V{1} = -3*phi{2} - 5*phi{3} - 3*phi{4};
    V{2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4};
    V{3} = -1/3*phi{2} - phi{3} - phi{4};
    
elseif ( strcmpi(schemeName, 'abnorsett5') == 1 )
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,1} = expinteg.psiEval(1, C(1), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute V:
    V{1} = -(4*phi{2} + 26/3*phi{3} + 9*phi{4} + 4*phi{5});
    V{2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5};
    V{3} = -(4/3*phi{2} + 14/3*phi{3} + 7*phi{4} + 4*phi{5});
    V{4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    
elseif ( strcmpi(schemeName, 'abnorsett6') == 1 )
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,1} = expinteg.psiEval(1, C(1), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute U:
    V{1} = -(5*phi{2} + 77/6*phi{3} + 71/4*phi{4} + 14*phi{5} + 5*phi{6});
    V{2} = 5*phi{2} + 107/6*phi{3} + 59/2*phi{4} + 26*phi{5} + 10*phi{6};
    V{3} = -(10/3*phi{2} + 13*phi{3} + 49/2*phi{4} + 24*phi{5} + 10*phi{6});
    V{4} = 5/4*phi{2} + 61/12*phi{3} + 41/4*phi{4} + 11*phi{5} + 5*phi{6};
    V{5} = -(1/5*phi{2} + 5/6*phi{3} + 7/4*phi{4} + 2*phi{5} + phi{6});
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ETD RUNGE-KUTTA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'etdrk2') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;

    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{2,1} = phi{1};
    
    % Compute B:
    B{1} = phi{1} - phi{2};
    B{2} = phi{2};
    
elseif ( strcmpi(schemeName, 'etdrk4') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = psi{1,2};
    A{4,3} = 2*psi{1,2};
    
    % Compute B:
    B{2} = 2*phi{2} - 4*phi{3};
    B{3} = 2*phi{2} - 4*phi{3};
    B{4} = -phi{2} + 4*phi{3};
    
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
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = psi{1,2};
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{1,5} = psi{1,2};
    psi{1,6} = expinteg.psiEval(1, C(6), LR, N, dim, nVars);
    psi{1,7} = expinteg.psiEval(1, C(7), LR, N, dim, nVars);
    psi{1,8} = phi{1};
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{2,4} = expinteg.psiEval(2, C(4), LR, N, dim, nVars);
    psi{2,6} = expinteg.psiEval(2, C(6), LR, N, dim, nVars);
    psi{2,7} = expinteg.psiEval(2, C(7), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    psi{3,6} = expinteg.psiEval(3, C(6), LR, N, dim, nVars);
    psi{3,7} = expinteg.psiEval(3, C(7), LR, N, dim, nVars);
    psi{4,6} = expinteg.psiEval(4, C(6), LR, N, dim, nVars);
    psi{4,7} = expinteg.psiEval(4, C(7), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*psi{2,2};
    A{4,3} = 2*psi{2,4};
    A{5,3} = -2*psi{2,2} + 16*psi{3,2};
    A{5,4} = 8*psi{2,2} - 32*psi{3,2};
    A{6,4} = 8*psi{2,6} - 32*psi{3,6};
    A{7,4} = -(125/162)*A{6,4};
    A{6,5} = -2*psi{2,6} + 16*psi{3,6};
    A{7,5} = (125/1944)*A{6,4} - (4/3)*psi{2,7} + (40/3)*psi{3,7};
    Phi = (5/32)*A{6,4} - (25/28)*psi{2,6} + (81/175)*psi{2,7} - ...
        (162/25)*psi{3,7} + (150/7)*psi{4,6} + (972/35)*psi{4,7} + 6*phi{4};
    A{8,5} = -(16/3)*phi{2} + (208/3)*phi{3} - 40*Phi;
    A{7,6} = (3125/3888)*A{6,4} + (25/3)*psi{2,7} - (100/3)*psi{3,7};
    A{8,6} = (250/21)*phi{2} - (250/3)*phi{3} + (250/7)*Phi;
    A{8,7} = (27/14)*phi{2} - 27*phi{3} + (135/7)*Phi;
    
    % Compute B:
    B{6} = (125/14)*phi{2} - (625/14)*phi{3} + (1125/14)*phi{4};
    B{7} = -(27/14)*phi{2} + (162/7)*phi{3} - (405/7)*phi{4};
    B{8} = (1/2)*phi{2} - (13/2)*phi{3} + (45/2)*phi{4};
    
elseif ( strcmpi(schemeName, 'friedli') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*psi{2,2};
    A{4,2} = -26/25*phi{1} + 2/25*phi{2};
    A{4,3} = 26/25*phi{1} + 48/25*phi{2};
    
    % Compute B:
    B{3} = 4*phi{2} - 8*phi{3};
    B{4} = -phi{2} + 4*phi{3};
    
elseif ( strcmpi(schemeName, 'hochbruck-ostermann') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    C(5) = 1/2;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{1,5} = expinteg.psiEval(1, C(5), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4*psi{2,2};
    A{4,2} = phi{2};
    A{5,2} = 1/4*phi{2} - phi{3} + 2*psi{2,2} - 4*psi{3,2};
    A{4,3} = A{4,2};
    A{5,3} = A{5,2};
    A{5,4} = psi{2,2} - A{5,2};
    
    % Compute B:
    B{4} = -phi{2} + 4*phi{3};
    B{5} = 4*phi{2} - 8*phi{3};
    
elseif ( strcmpi(schemeName, 'krogstad') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4*psi{2,2};
    A{4,3} = 2*phi{2};
    
    % Compute B:
    B{2} = 2*phi{2} - 4*phi{3};
    B{3} = 2*phi{2} - 4*phi{3};
    B{4} = -phi{2} + 4*phi{3};
    
elseif ( strcmpi(schemeName, 'minchev') == 1 ) 
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 4/25*psi{1,2} + 24/25*psi{2,2};
    A{4,2} = 21/5*phi{2} - 108/5*phi{3};
    A{4,3} = 1/20*phi{1} - 33/10*phi{2} + 123/5*phi{3};
    
    % Compute B:
    B{2} = -1/10*phi{1} + 1/5*phi{2} - 4*phi{3} + 12*phi{4};
    B{3} = 1/30*phi{1} + 23/5*phi{2} - 8*phi{3} - 4*phi{4};
    B{4} = 1/30*phi{1} - 7/5*phi{2} + 6*phi{3} - 4*phi{4};
    
elseif ( strcmpi(schemeName, 'strehmel-weiner') == 1 ) 
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 2*psi{2,2};
    A{4,2} = -2*phi{2};
    A{4,3} = 4*phi{2};
    
    % Compute B:
    B{3} = 4*phi{2} - 8*phi{3};
    B{4} = -phi{2} + 4*phi{3};
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAWSON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'ablawson4') == 1 )
    
    
    % Compute the phi-functions:
    phi0 = expinteg.phiEval(0, LR, N, dim, nVars);
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,1} = expinteg.psiEval(1, C(1), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        phi0 = real(phi0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    e2z = phi0.*phi0;
    e3z = e2z.*phi0;
    e4z = e2z.*e2z;
    
    % Compute B:
    B{1} = 55/24*phi0;
    
    % Compute V:
    V{1} = -59/24*e2z;
    V{2} = 37/24*e3z;
    V{3} = -3/8*e4z;
    
elseif ( strcmpi(schemeName, 'lawson4') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi0 = expinteg.phiEval(0, LR, N, dim, nVars);
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi0 = real(phi0);
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{2,1} = 1/2*psi02;
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{1} = 1/6*phi0;
    B{2} = 1/3*psi02;
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
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/6;
    
elseif ( strcmpi(schemeName, 'genlawson42') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -psi{2,2};
    U{3,1} = -psi{2,2} + 1/4;
    U{4,1} = -phi{2} + 1/2*psi02;
    
    % Compute V:
    V{1} = -phi{2} + 1/3*psi02 + 1/6;
    
elseif ( strcmpi(schemeName, 'genlawson43') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -2*psi{2,2} - 2*psi{3,2};
    U{3,1} = -2*psi{2,2} - 2*psi{3,2} + 5/8;
    U{4,1} = -2*phi{2} - 2*phi{3} + 5/4*psi02;
    U{2,2} = 1/2*psi{2,2} + psi{3,2};
    U{3,2} = 1/2*psi{2,2} + psi{3,2} - 3/16;
    U{4,2} = 1/2*phi{2} + phi{3} - 3/8*psi02;
    
    % Compute V:
    V{1} = -2*phi{2} - 2*phi{3} + 5/6*psi02 + 1/2;
    V{2} = 1/2*phi{2} + phi{3} - 1/4*psi02 - 1/6;
    
elseif ( strcmpi(schemeName, 'genlawson44') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    psi{4,2} = expinteg.psiEval(4, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -3*psi{2,2} - 5*psi{3,2} - 3*psi{4,2};
    U{2,2} = 3/2*psi{2,2} + 4*psi{3,2} + 3*psi{4,2};
    U{2,3} = -1/3*psi{2,2} - 1*psi{3,2} - 1*psi{4,2};
    U{3,1} = U{2,1} + 35/32;
    U{3,2} = U{2,2} - 21/32;
    U{3,3} = U{2,3} + 5/32;
    U{4,1} = -3*phi{2} - 5*phi{3} - 3*phi{4} + 35/16*psi02;
    U{4,2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4} - 21/16*psi02;
    U{4,3} = -1/3*phi{2} - phi{3} - phi{4} + 5/16*psi02;
    
    % Compute V:
    V{1} = -3*phi{2} - 5*phi{3} - 3*phi{4} + 35/24*psi02 + 1;
    V{2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4} - 7/8*psi02 - 2/3;
    V{3} = -1/3*phi{2} - 1*phi{3} - 1*phi{4} + 5/24*psi02 + 1/6;
    
elseif ( strcmpi(schemeName, 'genlawson45') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    psi{4,2} = expinteg.psiEval(4, C(2), LR, N, dim, nVars);
    psi{5,2} = expinteg.psiEval(5, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/6;
    
    % Compute U:
    U{2,1} = -4*psi{2,2} - 26/3*psi{3,2} - 9*psi{4,2} - 4*psi{5,2};
    U{2,2} = 3*psi{2,2} + 19/2*psi{3,2} + 12*psi{4,2} + 6*psi{5,2};
    U{2,3} = -4/3*psi{2,2} - 14/3*psi{3,2} - 7*psi{4,2} - 4*psi{5,2};
    U{2,4} = 1/4*psi{2,2} + 88/96*psi{3,2} + 3/2*psi{4,2} + psi{5,2};
    U{3,1} = U{2,1} + 105/64;
    U{3,2} = U{2,2} - 189/128;
    U{3,3} = U{2,3} + 45/64;
    U{3,4} = U{2,4} - 35/256;
    U{4,1} = -4*phi{2} - 26/3*phi{3} - 9*phi{4} - 4*phi{5} + 105/32*psi02;
    U{4,2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5} - 189/64*psi02;
    U{4,3} = -4/3*phi{2} - 14/3*phi{3} - 7*phi{4} - 4*phi{5} + 45/32*psi02;
    U{4,4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5} - 35/128*psi02;
    
    % Compute V:
    V{1} = -4*phi{2} - 26/3*phi{3} - 9*phi{4} - 4*phi{5} + 35/16*psi02 + ...
        5/3;
    V{2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5} - 63/32*psi02 - ...
        5/3;
    V{3} = -4/3*phi{2} - 14/3*phi{3} - 7*phi{4} - 4*phi{5} + ...
        15/16*psi02 + 5/6;
    V{4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5} - ...
        35/192*psi02 - 1/6;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFIED GENERALIZED LAWSON:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'modgenlawson41') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = phi{2} - 1/3*psi02;
    
elseif ( strcmpi(schemeName, 'modgenlawson42') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} =  1/2*phi{2} + phi{3} - 1/4*psi02;
    
    % Compute U:
    U{2,1} = -psi{2,2};
    U{3,1} = -psi{2,2} + 1/4;
    U{4,1} = -phi{2} + 1/2*psi02;
    
    % Compute V:
    V{1} = -1/2*phi{2} + phi{3} + 1/12*psi02;
    
elseif ( strcmpi(schemeName, 'modgenlawson43') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/3*phi{2} + phi{3} + phi{4} - 5/24*psi02;
    
    % Compute U:
    U{2,1} = -2*psi{2,2} - 2*psi{3,2};
    U{3,1} = -2*psi{2,2} - 2*psi{3,2} + 5/8;
    U{4,1} = -2*phi{2} - 2*phi{3} + 5/4*psi02;
    U{2,2} = 1/2*psi{2,2} + psi{3,2};
    U{3,2} = 1/2*psi{2,2} + psi{3,2} - 3/16;
    U{4,2} = 1/2*phi{2} + phi{3} - 3/8*psi02;
    
    % Compute V:
    V{1} = -phi{2} + phi{3} + 3*phi{4} + 5/24*psi02;
    V{2} = 1/6*phi{2} - phi{4} - 1/24*psi02;
    
elseif ( strcmpi(schemeName, 'modgenlawson44') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    psi{4,2} = expinteg.psiEval(4, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5} - 35/192*psi02;
    
    % Compute U:
    U{2,1} = -3*psi{2,2} - 5*psi{3,2} - 3*psi{4,2};
    U{2,2} = 3/2*psi{2,2} + 4*psi{3,2} + 3*psi{4,2};
    U{2,3} = -1/3*psi{2,2} - 1*psi{3,2} - 1*psi{4,2};
    
    U{3,1} = U{2,1} + 35/32;
    U{3,2} = U{2,2} - 21/32;
    U{3,3} = U{2,3} + 5/32;
    
    U{4,1} = -3*phi{2} - 5*phi{3} - 3*phi{4} + 35/16*psi02;
    U{4,2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4} - 21/16*psi02;
    U{4,3} = -1/3*phi{2} - phi{3} - phi{4} + 5/16*psi02;
    
    % Compute V:
    V{1} = -3/2*phi{2} + 1/2*phi{3} + 6*phi{4} + 6*phi{5} + 35/96*psi02;
    V{2} = 1/2*phi{2} + 1/3*phi{3} - 3*phi{4} - 4*phi{5} - 7/48*psi02;
    V{3} = -1/12*phi{2} - 1/12*phi{3} + 1/2*phi{4} + phi{5} + 5/192*psi02;
    
elseif ( strcmpi(schemeName, 'modgenlawson45') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1/2;
    C(3) = 1/2;
    C(4) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi02 = expinteg.psiEval(0, C(2), LR, N, dim, nVars);
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    psi{1,4} = expinteg.psiEval(1, C(4), LR, N, dim, nVars);
    psi{2,2} = expinteg.psiEval(2, C(2), LR, N, dim, nVars);
    psi{3,2} = expinteg.psiEval(3, C(2), LR, N, dim, nVars);
    psi{4,2} = expinteg.psiEval(4, C(2), LR, N, dim, nVars);
    psi{5,2} = expinteg.psiEval(5, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi02 = real(psi02);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/2;
    A{4,3} = psi02;
    
    % Compute B:
    B{2} = 1/3*psi02;
    B{3} = B{2};
    B{4} = 12/59*phi{2} + 50/59*phi{3} + 105/59*phi{4} + 120/59*phi{5} - ...
        60/59*phi{6} - 157/944*psi02;
    
    % Compute U:
    U{2,1} = -4*psi{2,2} - 26/3*psi{3,2} - 9*psi{4,2} - 4*psi{5,2};
    U{2,2} = 3*psi{2,2} + 19/2*psi{3,2} + 12*psi{4,2} + 6*psi{5,2};
    U{2,3} = -4/3*psi{2,2} - 14/3*psi{3,2} - 7*psi{4,2} - 4*psi{5,2};
    U{2,4} = 1/4*psi{2,2} + 88/96*psi{3,2} + 3/2*psi{4,2} + psi{5,2};
    U{3,1} = U{2,1} + 105/64;
    U{3,2} = U{2,2} - 189/128;
    U{3,3} = U{2,3} + 45/64;
    U{3,4} = U{2,4} - 35/256;
    U{4,1} = -4*phi{2} - 26/3*phi{3} - 9*phi{4} - 4*phi{5} + 105/32*psi02;
    U{4,2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5} - 189/64*psi02;
    U{4,3} = -4/3*phi{2} - 14/3*phi{3} - 7*phi{4} - 4*phi{5} + 45/32*psi02;
    U{4,4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5} - 35/128*psi02;
    
    % Compute V:
    V{1} = -116/59*phi{2} - 34/177*phi{3} + 519/59*phi{4} + ...
        964/59*phi{5} - 600/59*phi{6} + 495/944*psi02;
    V{2} = 57/59*phi{2} + 121/118*phi{3} - 342/59*phi{4} - ...
        846/59*phi{5} + 600/59*phi{6} - 577/1888*psi02;
    V{3} = -56/177*phi{2} - 76/177*phi{3} + 112/59*phi{4} + ...
        364/59*phi{5} - 300/59*phi{6} + 25/236*psi02;
    V{4} = 11/236*phi{2} + 49/708*phi{3} - 33/118*phi{4} - ...
        61/59*phi{5} + 60/59*phi{6} - 181/11328*psi02;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREDICTOR-CORRECTOR:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( strcmpi(schemeName, 'pec423') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute B:
    B{2} = 1/3*phi{2} + phi{3} + phi{4};
    
    % Compute U:
    U{2,1} = -2*phi{2} - 2*phi{3};
    U{2,2} = 1/2*phi{2} + phi{3};
    
    % Compute V:
    V{1} = -phi{2} + phi{3} + 3*phi{4};
    V{2} = 1/6*phi{2} - phi{4};
    
elseif ( strcmpi(schemeName, 'pecec433') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/3*phi{2} + phi{3} + phi{4};
    
    % Compute B:
    B{3} = 1/3*phi{2} + phi{3} + phi{4};
    
    % Compute U:
    U{2,1} = -2*phi{2} - 2*phi{3};
    U{3,1} = -phi{2} + phi{3} + 3*phi{4};
    U{2,2} = 1/2*phi{2} + phi{3};
    U{3,2} = 1/6*phi{2} - phi{4};
    
    % Compute V:
    V{1} = -phi{2} + phi{3} + 3*phi{4};
    V{2} = 1/6*phi{2} - phi{4};
    
elseif ( strcmpi(schemeName, 'pec524') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute B:
    B{2} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    
    % Compute U:
    U{2,1} = -3*phi{2} - 5*phi{3} - 3*phi{4};
    U{2,2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4};
    U{2,3} = -1/3*phi{2} - phi{3} - phi{4};
    
    % Compute V:
    V{1} = -3/2*phi{2} + 1/2*phi{3} + 6*phi{4} + 6*phi{5};
    V{2} = 1/2*phi{2} + 1/3*phi{3} - 3*phi{4} - 4*phi{5};
    V{3} = -1/12*phi{2} - 1/12*phi{3} + 1/2*phi{4} + phi{5};
    
elseif ( strcmpi(schemeName, 'pecec534') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    
    % Compute B:
    B{3} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    
    % Compute U:
    U{2,1} = -3*phi{2} - 5*phi{3} - 3*phi{4};
    U{2,2} = 3/2*phi{2} + 4*phi{3} + 3*phi{4};
    U{2,3} = -1/3*phi{2} - phi{3} - phi{4};
    
    U{3,1} = -3/2*phi{2} + 1/2*phi{3} + 6*phi{4} + 6*phi{5};
    U{3,2} = 1/2*phi{2} + 1/3*phi{3} - 3*phi{4} - 4*phi{5};
    U{3,3} = -1/12*phi{2} - 1/12*phi{3} + 1/2*phi{4} + phi{5};
    
    % Compute V:
    V{1} = -3/2*phi{2} + 1/2*phi{3} + 6*phi{4} + 6*phi{5};
    V{2} = 1/2*phi{2} + 1/3*phi{3} - 3*phi{4} - 4*phi{5};
    V{3} = -1/12*phi{2} - 1/12*phi{3} + 1/2*phi{4} + phi{5};
    
elseif ( strcmpi(schemeName, 'pec625') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute B:
    B{2} = 1/5*phi{2} + 5/6*phi{3} + 7/4*phi{4} + 2*phi{5} + phi{6};
    
    % Compute U:
    U{2,1} = -4*phi{2} - 26/3*phi{3} - 9*phi{4} - 4*phi{5};
    U{2,2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5};
    U{2,3} = -4/3*phi{2} - 14/3*phi{3} - 7*phi{4} - 4*phi{5};
    U{2,4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    
    % Compute V:
    V{1} = -2*phi{2} - 1/3*phi{3} + 17/2*phi{4} + 16*phi{5} + 10*phi{6};
    V{2} = phi{2} + 7/6*phi{3} - 11/2*phi{4} - 14*phi{5} - 10*phi{6};
    V{3} = -1/3*phi{2} - 1/2*phi{3} + 7/4*phi{4} + 6*phi{5} + 5*phi{6};
    V{4} = 1/20*phi{2} + 1/12*phi{3} - 1/4*phi{4} - phi{5} - phi{6};
    
elseif ( strcmpi(schemeName, 'pecec635') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/5*phi{2} + 5/6*phi{3} + 7/4*phi{4} + 2*phi{5} + phi{6};
    
    % Compute B:
    B{3} = 1/5*phi{2} + 5/6*phi{3} + 7/4*phi{4} + 2*phi{5} + phi{6};
    
    % Compute U:
    U{2,1} = -4*phi{2} - 26/3*phi{3} - 9*phi{4} - 4*phi{5};
    U{2,2} = 3*phi{2} + 19/2*phi{3} + 12*phi{4} + 6*phi{5};
    U{2,3} = -4/3*phi{2} - 14/3*phi{3} - 7*phi{4} - 4*phi{5};
    U{2,4} = 1/4*phi{2} + 11/12*phi{3} + 3/2*phi{4} + phi{5};
    U{3,1} = -2*phi{2} - 1/3*phi{3} + 17/2*phi{4} + 16*phi{5} + 10*phi{6};
    U{3,2} = phi{2} + 7/6*phi{3} - 11/2*phi{4} - 14*phi{5} - 10*phi{6};
    U{3,3} = -1/3*phi{2} - 1/2*phi{3} + 7/4*phi{4} + 6*phi{5} + 5*phi{6};
    U{3,4} = 1/20*phi{2} + 1/12*phi{3} - 1/4*phi{4} - phi{5} - phi{6};
    
    % Compute V:
    V{1} = -2*phi{2} - 1/3*phi{3} + 17/2*phi{4} + 16*phi{5} + 10*phi{6};
    V{2} = phi{2} + 7/6*phi{3} - 11/2*phi{4} - 14*phi{5} - 10*phi{6};
    V{3} = -1/3*phi{2} - 1/2*phi{3} + 7/4*phi{4} + 6*phi{5} + 5*phi{6};
    V{4} = 1/20*phi{2} + 1/12*phi{3} - 1/4*phi{4} - phi{5} - phi{6};
    
elseif ( strcmpi(schemeName, 'pec726') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    phi{7} = expinteg.phiEval(7, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute B:
    B{2} = 1/6*phi{2} + 137/180*phi{3} + 15/8*phi{4} + 17/6*phi{5} + ...
        5/2*phi{6} + phi{7};
    
    % Compute U:
    U{2,1} = -5*phi{2} - 77/6*phi{3} - 71/4*phi{4} - 14*phi{5} - 5*phi{6};
    U{2,2} = 5*phi{2} + 107/6*phi{3} + 59/2*phi{4} + 26*phi{5} + 10*phi{6};
    U{2,3} = -10/3*phi{2} - 13*phi{3} - 49/2*phi{4} - 24*phi{5} - 10*phi{6};
    U{2,4} = 5/4*phi{2} + 61/12*phi{3} + 41/4*phi{4} + 11*phi{5} + 5*phi{6};
    U{2,5} = -1/5*phi{2} - 5/6*phi{3} - 7/4*phi{4} - 2*phi{5} - phi{6};
    
    % Compute V:
    V{1} = -5/2*phi{2} - 17/12*phi{3} + 83/8*phi{4} + 57/2*phi{5} + ...
        65/2*phi{6} + 15*phi{7};
    V{2} = 5/3*phi{2} + 47/18*phi{3} - 8*phi{4} - 92/3*phi{5} - ...
        40*phi{6} - 20*phi{7};
    V{3} = -5/6*phi{2} - 19/12*phi{3} + 29/8*phi{4} + 37/2*phi{5} + ...
        55/2*phi{6} + 15*phi{7};
    V{4} = 1/4*phi{2} + 31/60*phi{3} - phi{4} - 6*phi{5} - 10*phi{6} - ...
        6*phi{7};
    V{5} = -1/30*phi{2} - 13/180*phi{3} + 1/8*phi{4} + 5/6*phi{5} + ...
        3/2*phi{6} + phi{7};
    
elseif ( strcmpi(schemeName, 'pecec736') == 1 )
    
    % Compute C:
    C(1) = 0;
    C(2) = 1;
    C(3) = 1;
    
    % Compute the phi-functions:
    phi{1} = expinteg.phiEval(1, LR, N, dim, nVars);
    phi{2} = expinteg.phiEval(2, LR, N, dim, nVars);
    phi{3} = expinteg.phiEval(3, LR, N, dim, nVars);
    phi{4} = expinteg.phiEval(4, LR, N, dim, nVars);
    phi{5} = expinteg.phiEval(5, LR, N, dim, nVars);
    phi{6} = expinteg.phiEval(6, LR, N, dim, nVars);
    phi{7} = expinteg.phiEval(7, LR, N, dim, nVars);
    
    % Compute the psi-functions:
    psi{1,2} = expinteg.psiEval(1, C(2), LR, N, dim, nVars);
    psi{1,3} = expinteg.psiEval(1, C(3), LR, N, dim, nVars);
    
    % Take real part for diffusive problems (real eigenvalues):
    if ( isreal(L) == 1 )
        phi = cellfun(@(f) real(f), phi, 'UniformOutput', 0);
        psi = cellfun(@(f) real(f), psi, 'UniformOutput', 0);
    end
    
    % Compute A:
    A{3,2} = 1/6*phi{2} + 137/180*phi{3} + 15/8*phi{4} + 17/6*phi{5} + ...
        5/2*phi{6} + phi{7};
    
    % Compute B:
    B{3} = 1/6*phi{2} + 137/180*phi{3} + 15/8*phi{4} + 17/6*phi{5} + ...
        5/2*phi{6} + phi{7};
    
    % Compute U:
    U{2,1} = -5*phi{2} - 77/6*phi{3} - 71/4*phi{4} - 14*phi{5} - 5*phi{6};
    U{2,2} = 5*phi{2} + 107/6*phi{3} + 59/2*phi{4} + 26*phi{5} + 10*phi{6};
    U{2,3} = -10/3*phi{2} - 13*phi{3} - 49/2*phi{4} - 24*phi{5} - 10*phi{6};
    U{2,4} = 5/4*phi{2} + 61/12*phi{3} + 41/4*phi{4} + 11*phi{5} + 5*phi{6};
    U{2,5} = -1/5*phi{2} - 5/6*phi{3} - 7/4*phi{4} - 2*phi{5} - phi{6};
    U{3,1} = -5/2*phi{2} - 17/12*phi{3} + 83/8*phi{4} + 57/2*phi{5} + ...
        65/2*phi{6} + 15*phi{7};
    U{3,2} = 5/3*phi{2} + 47/18*phi{3} - 8*phi{4} - 92/3*phi{5} - ...
        40*phi{6} - 20*phi{7};
    U{3,3} = -5/6*phi{2} - 19/12*phi{3} + 29/8*phi{4} + 37/2*phi{5} + ...
        55/2*phi{6} + 15*phi{7};
    U{3,4} = 1/4*phi{2} + 31/60*phi{3} - phi{4} - 6*phi{5} - 10*phi{6} - ...
        6*phi{7};
    U{3,5} = -1/30*phi{2} - 13/180*phi{3} + 1/8*phi{4} + 5/6*phi{5} + ...
        3/2*phi{6} + phi{7};
    
    % Compute V:
    V{1} = -5/2*phi{2} - 17/12*phi{3} + 83/8*phi{4} + 57/2*phi{5} + ...
        65/2*phi{6} + 15*phi{7};
    V{2} = 5/3*phi{2} + 47/18*phi{3} - 8*phi{4} - 92/3*phi{5} - ...
        40*phi{6} - 20*phi{7};
    V{3} = -5/6*phi{2} - 19/12*phi{3} + 29/8*phi{4} + 37/2*phi{5} + ...
        55/2*phi{6} + 15*phi{7};
    V{4} = 1/4*phi{2} + 31/60*phi{3} - phi{4} - 6*phi{5} - 10*phi{6} - ...
        6*phi{7};
    V{5} = -1/30*phi{2} - 13/180*phi{3} + 1/8*phi{4} + 5/6*phi{5} + ...
        3/2*phi{6} + phi{7};
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Put everything in SCHEMECOEFFS:
schemeCoeffs.A = A;
schemeCoeffs.B = B;
schemeCoeffs.C = C;
schemeCoeffs.U = U;
schemeCoeffs.V = V;

% Compute the missing oefficients using the summation properties of the coeffs:
% (Unless the scheme is ABLAWSON4, LAWSON4 or ETDRK2.)
if ( ~strcmpi(schemeName, 'lawson4') && ~strcmpi(schemeName, 'ablawson4') && ...
     ~strcmpi(schemeName, 'etdrk2') ) 
    schemeCoeffs = computeMissingCoeffs(K, schemeCoeffs, phi, psi);
end

% Precompute the E quantities:
E = cell(s+1, 1);
for i = 1:s
    E{i} = exp(C(i)*dt*L);
end
E{s+1} = exp(dt*L);
schemeCoeffs.E = E;

% Multiply by time-step:
schemeCoeffs.A = cellfun(@(A) dt*A, schemeCoeffs.A, 'UniformOutput', 0);
schemeCoeffs.B = cellfun(@(B) dt*B, schemeCoeffs.B, 'UniformOutput', 0);
schemeCoeffs.U = cellfun(@(U) dt*U, schemeCoeffs.U, 'UniformOutput', 0);
schemeCoeffs.V = cellfun(@(V) dt*V, schemeCoeffs.V, 'UniformOutput', 0);

end

function schemeCoeffs = computeMissingCoeffs(K, schemeCoeffs, phi, psi)
%COMPUTEMISSINGCOEFFS   Compute the missing oefficients of an EXPINTEG using
%the summation properties of the coefficients.
%   SCHEMECOEFFS = COMPUTEMISSINGCOEFFS(K, SCHEMECOEFFS, PHI, PSI) uses the row
%   summation properties to compute the SCHEMECOEFFS A{i,1} and B{1} of the
%   EXPINTEG K, using and the phi- and psi-functions.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the coefficients:
A = schemeCoeffs.A;
B = schemeCoeffs.B;
U = schemeCoeffs.U;
V = schemeCoeffs.V;

% Number of internal stages S and number of steps used Q:
s = K.stages;
q = K.steps;

% Compute the coefficients A{i,1} using the row summation property:
for i = 2:s
    A{i,1} = psi{1,i};
    for j = 2:i-1
        if ( isempty(A{i,j}) == 0 )
            A{i,1} = A{i,1} - A{i,j};
        end
    end
    for j = 1:q-1
        if ( isempty(U{i,j}) == 0 )
            A{i,1} = A{i,1} - U{i,j};
        end
    end
end

% Compute the coefficient B{1} using the row summation property:
B{1} = phi{1};
for i = 2:s
    if ( isempty(B{i}) == 0 )
        B{1} = B{1} - B{i};
    end
end
for i = 1:q-1
    if ( isempty(V{i}) == 0 )
        B{1} = B{1} - V{i};
    end
end

% Put everything in SCHEMECOEFFS:
schemeCoeffs.A = A;
schemeCoeffs.B = B;

end
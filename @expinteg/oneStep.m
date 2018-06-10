function [uSol, NuSol] = oneStep(K, ~, schemeCoeffs, Nc, Nv, nVars, S, uSol, NuSol)
%ONESTEP   Compute solution at t_{n+1} from solution at previous time-steps.
%   [USOL, NUSOL] = ONESTEP(K, SCHEMECOEFFS, NC, NV, NVARS, USOL, NUSOL) updates 
%   the solution USOL and its nonlinear evaluations NUSOL using the coefficients 
%   SCHEMECOEFFS of the EXPINTEG K, the nonlinear part of the PDE in coefficient
%   and value spaces NC and NV, the number of variables NVARS and the 
%   SPINOPERATOR S.

% Notes: USOL and NUSOL are cell-arrays of Fourier coefficients that represent 
% the solution at different time-steps. They have Q entries for a multi-step 
% scheme with Q steps. USOL{1} is the vector of Fourier coefficients of the 
% solution at t_{n} (a matrix in 2D and a tensor in 3D), USOL{2} is the solution 
% at t_{n-1}, ..., and USOL{q} is the solution at time t_{n-q+1}. The solution 
% at t_{n+1} is computed from the solution at previous time-steps and is stored 
% in USOL{1}, and the others are shifted from one position to the right, i.e.,
% USOL{2}=USOL{1}, USOL{3}=USOL{2}, ..., USOL{Q}=USOL{Q-1}. Same for NUSOL.
% The formula for an expnonential general linear method is:
%    
%   Internal stages: 
%
%   (1)     v_{1} = u_{n} and, for 2<=i<=s,
%           v_{i} = exp(c_{i}*h*L)*v_{1} + h*sum_{j=1}^{i-1} A_{i,j}*N(v_{j})
%                                        + h*sum_{j=1}^{q-1} U_{i,j}*N(u_{n-j})
%
%   Solution at t_{n+1}:
%
%   (2)     u_{n+1} = exp(h*L)*u_{n} + h*sum_{i=1}^{s} B_{i}*N(v_{i}) 
%                                    + h*sum_{i=1}^{q-1} V_{i}*N(u_{n-i})
%

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
s = K.stages;               % number of internal stages
q = K.steps;                % number of steps (>1 for multistep methods)
N = size(uSol{1}, 1)/nVars; % grid points

% Get the values to coeffs and coeffs to values transform:
v2c = getVals2CoeffsTransform(S);
c2v = getCoeffs2ValsTransform(S);

% Coefficients of the scheme:
A = schemeCoeffs.A;
B = schemeCoeffs.B;
E = schemeCoeffs.E;
U = schemeCoeffs.U;
V = schemeCoeffs.V;

% VSOL stores the internal stages and NVSOL stores the nonlinear evaluation of 
% these internal stages:
vSol = cell(s, 1);
NvSol = cell(s, 1);

% The first internal stage VSOL{1} is equal to USOL{1}, same for NVSOL:
vSol{1} = uSol{1};
NvSol{1} = NuSol{1};

% Compute the S-1 other internal stages:
for i = 2:s
    
    % Use formula (1):
    vSol{i} = E{i}.*vSol{1};
    for j = 1:i-1
        if ( isempty(A{i,j}) == 0 )
            vSol{i} = vSol{i} + A{i,j}.*NvSol{j};
        end
    end
    for j = 1:q-1
        if ( isempty(U{i,j}) == 0 )        
            vSol{i} = vSol{i} + U{i,j}.*NuSol{j+1};
        end
    end
    
    % Nonlinear evaluation of the internal stages with IFFT/FFT:
    vals = c2v(vSol{i}(1:N,:,:));
    for k = 1:nVars-1
        idx = k*N + 1;
        vals = [vals; c2v(vSol{i}(idx:idx+N-1,:,:))];
    end
    vals = Nv(vals);
    coeffs = v2c(vals(1:N,:,:));
    for k = 1:nVars-1
        idx = k*N + 1;
        coeffs = [coeffs; v2c(vals(idx:idx+N-1,:,:))];
    end
    NvSol{i} = Nc.*coeffs;
    
end

% Compute the solution at t_{n+1} with formula (2):
sol = E{s+1}.*vSol{1};
for i = 1:s
    if ( isempty(B{i}) == 0 )
        sol = sol + B{i}.*NvSol{i};
    end
end
for i = 1:q-1
    if ( isempty(V{i}) == 0 )
        sol = sol + V{i}.*NuSol{i+1};
    end
end

% Update the solution:
if ( q == 1 )
    uSol{1} = sol;
else
    uSol = {sol, uSol{1:end-1}}.';
end
 
% Nonlinear evaluation of the solution at t_{n+1}:
vals = c2v(uSol{1}(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    vals = [vals; c2v(uSol{1}(idx:idx+N-1,:,:))]; %#ok<*AGROW>
end
vals = Nv(vals);
coeffs = v2c(vals(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    coeffs = [coeffs; v2c(vals(idx:idx+N-1,:,:))];
end
Nsol = Nc.*coeffs;

% Update the nonlinear evaluations:
if ( q == 1 )
    NuSol{1} = Nsol;
else
    NuSol = {Nsol, NuSol{1:end-1}}.';
end

end
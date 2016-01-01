function uSol = oneStep(K, schemeCoeffs, Nc, S, uSol)
%ONESTEP   Compute solution at t_{n+1} from solution at previous timesteps.
%   USOL = ONESTEP(K, SCHEMECOEFFS, NC, S, USOL) updates the solution USOL using 
%   the coefficients SCHEMECOEFFS of the SPINSCHEME K, the nonlinear part of the 
%   PDE in coefficient space NC, and the SPINOPERATOR S.

% Notes: USOL is a cell-array of Fourier coefficients that represents the 
% solution at different timesteps. It has Q elements for a multi-step scheme 
% with Q steps. USOL{1} is the vector of Fourier coefficients of the solution at 
% t_{n} (a matrix in 2D and a tensor in 3D), USOL{2} is the solution at t_{n-1}, 
% ..., and USOL{q} is the solution at time t_{n-q+1}. The solution at t_{n+1}
% is computed from the solution at previous timesteps and is stored in USOL{1},
% and the others are shifted from one position to the right, i.e.,
% USOL{2}=USOL{1}, USOL{3}=USOL{2}, ..., USOL{Q}=USOL{Q-1}. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
A = schemeCoeffs.A;
B = schemeCoeffs.B;
E = schemeCoeffs.E;
Nv = S.nonlinearPartVals;
nVars = S.numVars;
q = K.steps;
s = K.internalStages;
U = schemeCoeffs.U;
V = schemeCoeffs.V;

% Number of grid points N:
N = size(uSol{1}, 1)/nVars;

% VSOL stores the internal stages and NL stores the nonlinear evaluation of 
% these internal stages:
vSol = cell(s, 1);
NvSol = cell(s, 1);

% The first internal stage VSOL{1} is equal to USOL{1}:
vSol{1} = uSol{1};

% Nonlinear evaluation of VSOL{1} in value space with IFFT/FFT:
vals = ifftn(vSol{1}(1:N,:,:));
for k = 1:nVars-1
   idx = k*N + 1;
   vals = [vals; ifftn(vSol{1}(idx:idx+N-1,:,:))]; %#ok<*AGROW>
end
vals = Nv(vals);
coeffs = fftn(vals(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
end
NvSol{1} = Nc.*coeffs;

% Compute the S-1 other internal stages:
for i = 2:s
    % v_{i} = exp(c_{i}*h*L)*v_{1} + ...
    vSol{i} = E{i}.*vSol{1};
    % ... + h*sum_{j=1}^{i-1}A_{i,j}*N(v_{j}) + ...
    for j = 1:i-1
        if ( isempty(A{i,j}) == 0 )
            vSol{i} = vSol{i} + A{i,j}.*NvSol{j};
        end
    end
    % ... + h*sum_{j=1}^{q-1}U_{i,j}*N(u_{n-j})
    for j = 1:q-1
        if ( isempty(U{i,j}) == 0 )        
            % Nonlinear evaluation of USOL{J+1}=u_{n-j} with IFFT/FFT:
            vals = ifftn(uSol{j+1}(1:N,:,:));
            for k = 1:nVars-1
                idx = k*N + 1;
                vals = [vals; ifftn(uSol{j+1}(idx:idx+N-1,:,:))];
            end
            vals = Nv(vals);
            coeffs = fftn(vals(1:N,:,:));
            for k = 1:nVars-1
                idx = k*N + 1;
                coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
            end
            vSol{i} = vSol{i} + U{i,j}.*(Nc.*coeffs);
        end
    end
    % Nonlinear evaluation of VSOL{I}=v_{i} with IFFT/FFT:
    vals = ifftn(vSol{i}(1:N,:,:));
    for k = 1:nVars-1
        idx = k*N + 1;
        vals = [vals; ifftn(vSol{i}(idx:idx+N-1,:,:))];
    end
    vals = Nv(vals);
    coeffs = fftn(vals(1:N,:,:));
    for k = 1:nVars-1
        idx = k*N + 1;
        coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
    end
    NvSol{i} = Nc.*coeffs;
end

% Compute the solution at t_{n+1}:
% u_{n+1} = exp(h*L)*u_{n} + ...
vSol{1} = E{s+1}.*vSol{1};
% ... + h*sum_{i=1}^{s}B_{i}*N(v_{i})
for i = 1:s
    if ( isempty(B{i}) == 0 )
        vSol{1} = vSol{1} + B{i}.*NvSol{i};
    end
end
% ... + h*sum_{i=1}^{q-1}V_{i}*N(u_{n-i})
for i = 1:q-1
    if ( isempty(V{i}) == 0 )
        % Nonlinear evaluation of USOL{I+1}=u_{n-i} with IFFT/FFT:
        vals = ifftn(uSol{i+1}(1:N,:,:));
        for k = 1:nVars-1
            idx = k*N + 1;
            vals = [vals; ifftn(uSol{i+1}(idx:idx+N-1,:,:))];
        end
        vals = Nv(vals);
        coeffs = fftn(vals(1:N,:,:));
        for k = 1:nVars-1
            idx = k*N + 1;
            coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
        end
        vSol{1} = vSol{1} + V{i}.*(Nc.*coeffs);
    end
end

% Update the solution:
if ( q == 1 )
    uSol{1} = vSol{1};
else
    uSol = {vSol{1}, uSol{1:end-1}}.';
end
 
end
function [uSol, NuSol] = oneStep(K, dt, schemeCoeffs, gc, gv, nVars, S, uSol, NuSol)
%ONESTEP   Compute solution at t_{n+1} from solution at previous time-steps.
%   [USOL, NUSOL] = ONESTEP(K, SCHEMECOEFFS, GC, GV, NVARS, USOL, NUSOL) updates 
%   the solution USOL and its nonlinear evaluations NUSOL using the coefficients 
%   SCHEMECOEFFS of the IMEX K, the nonlinear parts of the PDE in coefficient 
%   and value space GC and GV, the number of variables NVARS, the SPINOPERATOR S 
%   and the time-step DT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note (1): USOL and NUSOL are cell-arrays of Fourier coeffs that represent 
% the solution at different time-steps. They have Q entries for a multi-step 
% scheme with Q steps. USOL{1} is the vector of Fourier coefficients of the 
% solution at t_{n}, USOL{2} is the solution at t_{n-1}, ..., and USOL{q} is 
% the solution at time t_{n-q+1}. The solution at t_{n+1} is computed from the 
% solution at previous time-steps and is stored in USOL{1}, and the others are 
% shifted from a position to the right, i.e., USOL{2}=USOL{1}, USOL{3}=USOL{2}, 
% ..., USOL{Q}=USOL{Q-1}. Same for NUSOL. 

% Note (2): For the moment, we only support the (one-step) LIRK4 scheme and 
% PDEs with linear part = Laplacian.

% Set-up:
q = K.steps;                      % number of steps (>1 for multistep methods)
N = sqrt(size(uSol{1}, 1)/nVars); % number of grid points

% Get the values to coeffs and coeffs to values transform:
v2c = getVals2CoeffsTransform(S);
c2v = getCoeffs2ValsTransform(S);

% The values and the nonlinear evaluations are stored in USOL and NUSOL:
v = uSol{1};
Nv = NuSol{1};

% Extract the coefficients of the scheme:
P = schemeCoeffs.precond;
Lap = schemeCoeffs.linmat;

% One step of LIRK4:
if ( strcmpi(K.scheme, 'lirk4') == 1 )
    L = schemeCoeffs.lufactors{1, 1};
    U = schemeCoeffs.lufactors{2, 1};
    La = schemeCoeffs.lufactors{1, 2};
    Ua = schemeCoeffs.lufactors{2, 2};
    w = P*v;
    wa = w + dt*P*1/4*Nv;
    a = Ua\(La\wa);
    Na = gc*vals2coeffs(gv(coeffs2vals(a, N, nVars, c2v)), N, nVars, v2c);
    wb = w + dt*Lap*1/2*a + dt*P*(-1/4*Nv + Na);
    b = Ua\(La\wb);
    Nb = gc*vals2coeffs(gv(coeffs2vals(b, N, nVars, c2v)), N, nVars, v2c);
    wc = w + dt*Lap*(17/50*a - 1/25*b) + dt*P*(-13/100*Nv + 43/75*Na + 8/75*Nb);
    c = Ua\(La\wc);
    Nc = gc*vals2coeffs(gv(coeffs2vals(c, N, nVars, c2v)), N, nVars, v2c);
    wd = w + dt*Lap*(371/1360*a - 137/2720*b + 15/544*c) ...
        + dt*P*(-6/85*Nv + 42/85*Na + 179/1360*Nb - 15/272*Nc);
    d = Ua\(La\wd);
    Nd = gc*vals2coeffs(gv(coeffs2vals(d, N, nVars, c2v)), N, nVars, v2c);
    we = w + dt*Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d) ...
        + dt*P*(79/24*Na - 5/8*Nb + 25/2*Nc - 85/6*Nd);
    e = Ua\(La\we);
    Ne = gc*vals2coeffs(gv(coeffs2vals(e, N, nVars, c2v)), N, nVars, v2c);
    v = v + dt*(U\(L\(Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d + 1/4*e))))...
        + dt*(25/24*Na - 49/48*Nb + 125/16*Nc - 85/12*Nd + 1/4*Ne);

% One step of IMEXBDF4:
elseif ( strcmpi(K.scheme, 'imexbdf4') == 1 )
    L = schemeCoeffs.lufactors{1, 1};
    U = schemeCoeffs.lufactors{2, 1};
    v = U\(L\(P*(48*uSol{1} - 36*uSol{2} + 16*uSol{3} - 3*uSol{4} + ...
        + 48*dt*NuSol{1} - 72*dt*NuSol{2} + 48*dt*NuSol{3} - 12*dt*NuSol{4})));
end

% Nonlinear evaluation of the solution at t_{n+1}:
Nv = gc*vals2coeffs(gv(coeffs2vals(v, N, nVars, c2v)), N, nVars, v2c);

% Update the solutions USOL:
if ( q == 1 )
    uSol{1} = v;
else
    uSol = {v, uSol{1:end-1}}.';
end

% Update the nonlinear evaluations NUSOL:
if ( q == 1 )
    NuSol{1} = Nv;
else
    NuSol = {Nv, NuSol{1:end-1}}.';
end

end

function vals = coeffs2vals(coeffs, N, nVars, c2v)
%COEFFS2VALS   Coeffs to values transform with NVARS variables.

vals = c2v(coeffs(1:N^2));
for k = 1:nVars-1
    idx = k*N^2 + 1;
    vals = [vals; c2v(coeffs(idx:idx+N^2-1))]; %#ok<*AGROW>
end

end

function coeffs = vals2coeffs(vals, N, nVars, v2c)
%VALS2COEFFS  Values to coeffs transform with NVARS variables.

coeffs = v2c(vals(1:N,:)); 
for k = 1:nVars-1
    idx = k*N + 1;
    coeffs = [coeffs; v2c(vals(idx:idx+N-1,:))];
end

end
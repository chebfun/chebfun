function [uSol, NuSol] = oneStep(~, dt, schemeCoeffs, gc, gv, ~, S, uSol, NuSol)
%ONESTEP   Compute solution at t_{n+1} from solution at previous time-steps.
%   [USOL, NUSOL] = ONESTEP(K, SCHEMECOEFFS, GC, GV, NVARS, USOL, NUSOL) updates 
%   the solution USOL and its nonlinear evaluations NUSOL using the coefficients 
%   SCHEMECOEFFS of the SPINSCHEME K, the nonlinear parts of the PDE in 
%   coefficient and value space GC and GV, the number of variables NVARS, the
%   SPINOPERATOR S and the time-step DT.

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

% Note (2): For the moment, we only support the LIRK4 scheme.

% Get the values to coeffs and coeffs to values transform:
v2c = getVals2CoeffsTransform(S);
c2v = getCoeffs2ValsTransform(S);

% Coefficients of the scheme:
Tsin2 = schemeCoeffs.multmat;
Lap = schemeCoeffs.L;
L = schemeCoeffs.LU{1, 1};
U = schemeCoeffs.LU{2, 1};
La = schemeCoeffs.LU{1, 2};
Ua = schemeCoeffs.LU{2, 2};
  
v = uSol{1};
Nv = NuSol{1};

w = Tsin2*v;
wa = w + dt*Tsin2*1/4*Nv;
a = Ua\(La\wa); 
Na = gc*v2c(gv(c2v(a)));
wb = w + dt*Lap*1/2*a + dt*Tsin2*(-1/4*Nv + Na);
b = Ua\(La\wb); 
Nb = gc*v2c(gv(c2v(b)));
wc = w + dt*Lap*(17/50*a - 1/25*b) + dt*Tsin2*(-13/100*Nv + 43/75*Na + 8/75*Nb);
c = Ua\(La\wc); 
Nc = gc*v2c(gv(c2v(c)));
wd = w + dt*Lap*(371/1360*a - 137/2720*b + 15/544*c) ...
    + dt*Tsin2*(-6/85*Nv + 42/85*Na + 179/1360*Nb - 15/272*Nc);
d = Ua\(La\wd); 
Nd = gc*v2c(gv(c2v(d)));
we = w + dt*Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d) ...
    + dt*Tsin2*(79/24*Na - 5/8*Nb + 25/2*Nc - 85/6*Nd);
e = Ua\(La\we); 
Ne = gc*v2c(gv(c2v(e)));
v = v + dt*(U\(L\(Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d + 1/4*e)))) ...
    + dt*(25/24*Na - 49/48*Nb + 125/16*Nc - 85/12*Nd + 1/4*Ne);

uSol{1} = v;
NuSol{1} = gc*v2c(gv(c2v(v)));

end
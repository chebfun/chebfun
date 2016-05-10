function [uSol, NuSol, dt] = startMultistep(K, adaptiveTime, dt, L, Nc, Nv, ...
    pref, S, uSol, NuSol)
%STARTMULTISTEP  Get enough initial data when using a multistep scheme.
%    [USOL, NUSOL, DT] = STARTMULTISTEP(K, ADAPTIVETIME, DT, L, NC, NV, ...
%    PREF, S, USOL, NUSOL) does a few steps of a one-step scheme with time-step 
%    DT to get enough initial data to start the multistep SPINSCHEME K, using 
%    the linear part L, the nonlinear parts of the operator in coefficient and
%    value space NC and NV, the SPINPREFERENCE object PREF, and the SPINOPERATOR
%    S. ADAPTIVETIME is 1 if adpative in time, 0 otherwise.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
errTol = pref.errTol;       % error tolerance 
M = pref.M;                 % points for the contour integral
q = K.steps;                % number of steps 
nVars = S.numVars;          % number of unknown functions

% Create a cell-array to store the coefficients at the Q steps:
coeffs = cell(q, 1);
Ncoeffs = cell(q, 1);

% Store the initial conidition in the last column:
coeffs{q} = uSol{1};
Ncoeffs{q} = NuSol{1};

% Set-up the scheme:
K = spinscheme('etdrk4');
schemeCoeffs = computeCoeffs(K, dt, L, M, S);
if ( adaptiveTime == 1 )
    schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
end

% Do Q-1 steps:
iter = 1;
while ( iter <= q-1 ) 
    
    [cNew, NcNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uSol, NuSol);
     
    % If adpative in time, two steps in time with DT/2:
    if ( adaptiveTime == 1 )
        [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, uSol, NuSol);
        [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, cNew2, ...
            NcNew2);
        err = max(abs(cNew{1}(:) - cNew2{1}(:)));
        err = err/max(abs(cNew2{1}(:)));
            
        % If successive step, store it:
        if ( err < errTol ) 
            coeffs{q-iter} = cNew2{1};
            Ncoeffs{q-iter} = NcNew2{1};
            uSol = cNew2;
            iter = iter + 1;
            
        % Otherwise, redo all the steps with DT/2:
        else
            dt = dt/2;
            schemeCoeffs = schemeCoeffs2;
            schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
            uSol = [];
            NuSol = [];
            uSol{1} = coeffs{q};
            NuSol{1} = Ncoeffs{q};
            iter = 1;
        end
        
    % If not adaptive in time, keep CNEW:
    else
        coeffs{q-iter} = cNew{1};
        Ncoeffs{q-iter} = NcNew{1};
        uSol = cNew;
        NuSol = NcNew;
        iter = iter + 1;
    end

end
uSol = coeffs;
NuSol = Ncoeffs;

end
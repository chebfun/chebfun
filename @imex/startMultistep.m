function [uSol, NuSol] = startMultistep(K, dt, L, Nc, Nv, ~, S, uInit, NuInit)
%STARTMULTISTEP   Get enough initial data when using a multistep scheme.
%    [USOL, NUSOL] = STARTMULTISTEP(K, dt, L, NC, NV, pref, S, uInit, NuInit)
%    uses a one-step algorithm to start a multistep scheme.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers: multistep algorithms are started with the one-step LIRK4
% time-stepping scheme.

% Set-up:
q = K.steps;                      % number of steps 
nVars = S.numVars;                % number of unknown functions
%N = sqrt(size(uSol{1}, 1)/nVars); % number of grid points

% Create a cell-array to store the coefficients at the Q steps:
uSol = cell(q, 1);
NuSol = cell(q, 1);

% Store the initial conidition in the last column:
uSol{q} = uInit{1};
NuSol{q} = NuInit{1};

% Set-up the scheme:
K = imex('lirk4');
schemeCoeffs = computeCoeffs(K, dt, L, [], S);

% Do (Q-1) steps of LIRK4:
uOld = uInit;
NuOld = NuInit;
for j = 1:(q-1)
    [uNew, NuNew] = oneStep(K, dt, schemeCoeffs, Nc, Nv, nVars, S, uOld, NuOld);
    uSol{q-j} = uNew{1};
    NuSol{q-j} = NuNew{1};
    uOld = uNew;
    NuOld = NuNew;
end

end
function [uSol, NuSol] = startMultistep(K, dt, L, Nc, Nv, pref, S, uInit, NuInit)
%STARTMULTISTEP  Get enough initial data when using a multistep scheme.
%    [USOL, NUSOL] = STARTMULTISTEP(K, dt, L, NC, NV, pref, S, uInit, NuInit)
%    uses a one-step algorithm with time-step DT, combined with a fixed point
%    algorithm, to get enough initial data to start the multistep SPINSCHEME K
%    using the linear part L, the nonlinear parts of the operator in coeff and
%    value space NC and NV, the SPINPREFERENCE object PREF, and the SPINOPERATOR
%    S.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers. The algorithm is in two parts:
% 1. Get an approximation of the data with ETDRK2.
% 2. Use a fixed point algorithm to refine this approximation.

%% Part 1: ETDRK2.

% Set-up:
M = pref.M;                 % points for the contour integral
q = K.steps;                % number of steps 
%errTol = max(1e-12, dt^q);  % error tolerance   
errTol = 1e-10;
nVars = S.numVars;          % number of unknown functions
N = size(L, 1)/nVars;       % grid points
dim = getDimension(S);      % spatial dimension (1, 2 or 3)

% Create a cell-array to store the coefficients at the Q steps:
uSol = cell(q, 1);
NuSol = cell(q, 1);

% Store the initial conidition in the last column:
uSol{q} = uInit{1};
NuSol{q} = NuInit{1};

% Set-up the scheme:
K = spinscheme('exprk5s8');
h = 10;
schemeCoeffs = computeCoeffs(K, dt/h, L, M, S);

% Create a contour around each eigenvalue of the linear part L:
LR = computeLR(S, dt, L, M);
%r = exp(2i*pi*((1:M) - .5)/M);
%LR = dt*repmat(L(:), 1, M) + repmat(r, nVars*N^dim, 1);

% Do Q-1 steps:
% uOld{1} = uInit{1};
% NuOld{1} = NuInit{1};
% for j = 1:q-1
%     [uNew, NuNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uOld, NuOld);
%     uSol{q-j} = uNew{1};
%     NuSol{q-j} = NuNew{1};
%     uOld = uNew;
%     NuOld = NuNew;
% end
% uOld{1} = uInit{1};
% NuOld{1} = NuInit{1};
% for k = 1:h*(q-1)
%     [uNew, NuNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uOld, NuOld);
%     uOld{1} = uNew{1};
%     NuOld{1} = NuNew{1};
%     if ( mod (k, h) == 0 )
%         uSol{q-k/h} = uNew{1};
%         NuSol{q-k/h} = NuNew{1};
%     end
% end
%uSol{1} = uNew{1};
%NuSol{1} = NuNew{1};

% uOld{1} = uInit{1};
% NuOld{1} = NuInit{1};
% [uNew, NuNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uOld, NuOld);
% uSol{1} = uNew{1};
% NuSol{1} = NuNew{1};

%% Part 2: Fixed point algorithm.

% Get the gamma-functions:
g = cell(q-1);
g0 = cell(q-1, 1);
for j = 1:q-1
    g0{j} = spinscheme.gammaEval(0, j, LR, N, dim, nVars);
    for k = 1:q-1
        g{j, k} = spinscheme.gammaEval(j, k, LR, N, dim, nVars);
    end
end
%phi1 = spinscheme.phiEval(1, LR, N, dim, nVars);
%phi2 = spinscheme.phiEval(2, LR, N, dim, nVars);
%g{1,1} = phi2;

% Take real part for diffusive problems (real eigenvalues):
if ( isreal(L) == 1 )
    g = cellfun(@(f) real(f), g, 'UniformOutput', 0);
    g0 = cellfun(@(f) real(f), g0, 'UniformOutput', 0);
end

% err = 1;
% iter = 0;
% iterMax = 100;
% uOld = uSol;
% NuOld = NuSol;
% while ( err > errTol && iter <= iterMax )
%     New coefficients:
%     for j = 1:q-1
%         ksi = exp(j*dt*L).*uOld{q} + dt*g0{j}.*NuOld{q};
%         uNew{q-j} = ksi;
%         for k = 1:q-1
%             uNew{q-j} = uNew{q-j} + dt*g{j, k}.*NuOld{q-k};
%         end
%     end
%     uNew{q} = uOld{q};
%     Nonlinear evaluations of the coefficients:
%     err = 0;
%     for j = 1:q-1
%         vals = ifftn(uNew{j}(1:N,:,:));
%         for k = 1:nVars-1
%             idx = k*N + 1;
%             vals = [vals; ifftn(uNew{j}(idx:idx+N-1,:,:))]; %#ok<*AGROW>
%         end
%         vals = Nv(vals);
%         coeffs = fftn(vals(1:N,:,:));
%         for k = 1:nVars-1
%             idx = k*N + 1;
%             coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
%         end
%         NuNew{j} = Nc.*coeffs;
%         err = max(err, norm(uOld{j} - uNew{j}, inf));
%     end
%     NuNew{q} = NuOld{q};
%     uOld = uNew;
%     NuOld = NuNew;
%     iter = iter + 1
%     err
% end
% uSol = uNew;
% NuSol = NuNew;

% err = 1;
% iter = 0;
% iterMax = 20;
% uOld{1} = uSol{1};
% NuOld{1} = NuSol{1};
% while ( err > errTol && iter < iterMax )
%     uNew{1} = exp(dt*L).*uInit{1} + dt*g0{1}.*NuInit{1} - dt*g{1,1}.*NuInit{1} + ...
%         + dt*g{1,1}.*NuOld{1};
%     %uNew{1} = uOld{q} + dt*phi1.*(L.*uOld{q} + NuOld{q}) - dt*g{1,1}.*NuOld{q} + ...
%     %    + dt*g{1,1}.*NuOld{1};
%     err = norm(uOld{1} - uNew{1}, inf);
%     NuNew{1} = Nc.*fft(Nv(ifft(uNew{1})));
%     uOld{1} = uNew{1};
%     NuOld{1} = NuNew{1};
% 	iter = iter + 1;
% end
% uSol{1} = uNew{1};
% NuSol{1} = NuNew{1};
% uSol{q} = uInit{1};
% NuSol{q} = NuInit{1};

uOld{1} = uInit{1};
NuOld{1} = NuInit{1};
iterMax = 20;
for k = 1:h*(q-1)
    [uNew, NuNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, uOld, NuOld);
    uOld{1} = uNew{1};
    NuOld{1} = NuNew{1};
    err = 1;
    iter = 0;
    while ( err > errTol && iter < iterMax )
        uNew{1} = exp(dt*L).*uInit{1} + dt*g0{1}.*NuInit{1} - ...
            dt*g{1,1}.*NuInit{1} + dt*g{1,1}.*NuOld{1};
        err = norm(uOld{1} - uNew{1}, inf);
        NuNew{1} = Nc.*fft(Nv(ifft(uNew{1})));
        uOld{1} = uNew{1};
        NuOld{1} = NuNew{1};
        iter = iter + 1;
    end
    if ( mod (k, h) == 0 )
        uSol{q-k/h} = uNew{1};
        NuSol{q-k/h} = NuNew{1};
    end
end

% options = optimoptions('fsolve');
% options.TolFun = errTol;
% options.TolX = errTol;
% u1 = fsolve(@(u1) myfun(u1, dt, L, uSol{q}, NuSol{q}, g0{1}, g{1,1}, Nc, Nv), ...
%     uSol{1}, options);
% uSol{1} = u1;
% NuSol{1} = Nc.*fft(Nv(ifft(u1)));

end

% function y = myfun(u1, dt, L, u0, Nu0, g01, g11, Nc, Nv)
%     y = u1 - exp(dt*L).*u0 - dt*g01.*Nu0 + dt*g11.*Nu0 - ...
%         dt*g11.*(Nc.*fft(Nv(ifft(u1))));
% end
function [u, info] = solvebvpNonlinear(N, rhs, L, u0, res, pref)
% Solve a nonlinear BVP, using damped Newton iteration in function space.
% Here
%   N:      Nonlinear chebop
%   rhs:    Right hand side of ODE
%   L:      Linearisation of N around the initial guess
%   u0:     Initial guess of solution
%   res:    Residual of ODE at initial guess
%   pref:   Cheboppref preference structure
%
% References: TODO: Fill out info for Deuflhard, DPhil thesis, and other
% references as appropriate


% Store preferences used in the Newton iteration in separate variables
maxIter  = pref.maxIter;
discType = pref.discretization;
delTol   = pref.deltol;
% Did the user request damped or undamped Newton iteration?
damped = pref.damped;
% Minimum allowed steplength
lambdaMin = pref.lambdaMin;
% NUMVARS indicate how many unknown function we seek.
numVars = nargin(N.op) - 1;

% Assign initial guess to u
u = u0;

% Store the domain we're working with.
dom = N.domain;

% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Print info to command window, and/or show plot of progress
[displayFig, displayTimer] = N.displayInfoInit(u0, pref);

% u = u + delta;
% ub = u.blocks;
% res = N.op(x, ub{:}) - rhs;

% Counter for number of Newton steps taken.
newtonCounter = 0;

% Variable that controls whether we want to stop the Newton iteration, either
% because we converged, or the process has been identified to be nonconvergent.
terminate = 0;

% Store a vector with information about the norm of the Newton updates
normDeltaVec = zeros(maxIter,1);

% Initial damping parameter
lambda = 1;
while ( ~terminate )
    % Compute a Newton update
    delta = -(L\res);
    
    % Store the norm of the update
    normDelta = mynorm(delta);
    % TODO: Old variable name, to be removed
    nrmDelta = normDelta;
    
    % Are we in damped mode?
    if ( damped )
        % Find a step-length lambda using an affine covariant algorithm due to
        % Deuflhard [!!!reference], also described in [!!!thesis]
        
        accept = 0;
        initPrediction = 1;
        while ~accept
            
            if newtonCounter > 0 && initPrediction
                mu = (nrmDeltaOld*nrmDeltaBar)/(mynorm(deltaBar-delta,'fro')*nrmDelta)*lambda;
                lambda = min(1,mu);
                initPrediction = 0;
            end
            
            if lambda < lambdaMin
                disp('Convergence failure')
                % Take full Newton step
                u = u + delta;
                accept = 1;
                nrmDeltaOld = nrmDelta;
                initPrediction = 1;
                lambda = 1;
                continue
            end
            
            uTrial = u + lambda*delta;
            
            uTrialb = uTrial.blocks;
            
            deResFunTrial = N.op(x, uTrialb{:}) - rhs;
            deltaBar = -(L\deResFunTrial);
            
            nrmDeltaBar = mynorm(deltaBar);
            
            cFactor = nrmDeltaBar/nrmDelta;
            
            muPrime = (.5*nrmDelta*lambda^2)/(mynorm(deltaBar-(1-lambda)*delta,'fro'));
            
            if cFactor >=1
                lambda = min(muPrime,.5*lambda);
                continue;
            end
            
            lambdaPrime = min(1,muPrime);
            
            if lambdaPrime == 1 && nrmDeltaBar < delTol
                u = uTrial + deltaBar;
                newtonCounter
                terminate = 1;
                break
            end
            
            % Switch to pure Newton if we are experiencing good convergence
            if lambdaPrime == 1 && cFactor < .5
                damped = 0;
            end
            
            
            if lambdaPrime >= 4*lambda
                lambda = lambdaPrime;
                continue;
            end
            
            % If we get all the way here, accept iterate
            accept = 1;
        end
        u = uTrial;
    else    % We are in undamped phase
        u = u + delta;
        cFactor = normDelta/normDeltaOld;
    end
    
    % Evaluate the residual
    ub = u.blocks;
    res = N.op(x, ub{:}) - rhs;
    
    % Update counter of Newton steps taken
    newtonCounter = newtonCounter +1;
        
    % Store information about the norm of the updates
    normDeltaVec(newtonCounter) = normDelta;
    nrmDeltaOld = nrmDelta; % Useful for damping strategy
    normDeltaOld = normDelta;
    % Print info to command window, and/or show plot of progress
    N.displayInfoIter(u, delta, newtonCounter, normDelta, cFactor, ...
        length(delta{1}), lambda, displayFig, displayTimer, pref)
    
    % TODO: Replace with error estimate -- introduce errorTol in cheboppref
    if ( normDelta < delTol )
        terminate = 1;
        %                     elseif (newt > 3 && normUpdate(newt) > 0.1*normUpdate(newt-3))
        %                         warning('CHEBFUN:bvpsc','Newton iteration stagnated.')
        %                         break
    elseif (newtonCounter > maxIter)
        warning('CHEBOP:solvebvpNonlinear','Newton iteration failed.')
        break
    else
        % Linearise around current solution:
        L = linearise(N, x, ub, []); % flag to negate contraint RHSs.
        
        % Assign the preferred discretisation method to the linop.
        L.discretization = discType;
    end
end

% Return useful information in the INFO structure
info.normDelta = normDeltaVec(1:newtonCounter);

end

function out = mynorm(f,type)
if ( isa(f, 'chebmatrix') )
    out = max(cellfun(@(u) get(u, 'vscale'), f.blocks));
else
    out = get(f, 'vscale');
end
end
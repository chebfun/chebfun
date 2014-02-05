function [u, info] = solvebvpNonlinear(N, rhs, L, u0, res, pref, displayInfo)
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

% TODO: Replace with real values further down in code (temp values to test
% display methods)
errEst = NaN;
errEstDE = NaN;
errEstBC = .5678;
% Store preferences used in the Newton iteration in separate variables
maxIter  = pref.maxIter;
discType = pref.discretization;
errTol   = pref.errTol;
% Did the user request damped or undamped Newton iteration?
damped = pref.damped;
% Minimum allowed steplength
lambdaMin = pref.lambdaMin;

% Assign initial guess to u
u = u0;

% Store the domain we're working with.
dom = N.domain;

% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Print info to command window, and/or show plot of progress
[displayFig, displayTimer] = displayInfo('init', u0, pref);

% Counter for number of Newton steps taken.
newtonCounter = 0;

% Variable that controls whether we want to stop the Newton iteration, either
% because we converged, or the process has been identified to be nonconvergent.
terminate = 0;

% Store a vector with information about the norm of the Newton updates
normDeltaVec = zeros(maxIter,1);

% Initial damping parameter
lambda = 1;

% Need to subtract the rhs from the residual passed in
res = res - rhs;
while ( ~terminate )
    % Compute a Newton update
    [delta, disc] = linsolve(L, res, discType);
    % We had two output arguments above, need to negate delta
    delta = -delta;
    %     delta = -(L\res);     % Old fashioned
    % Store the norm of the update
    normDelta = mynorm(delta);
    
    % Are we in damped mode?
    if ( damped )
        % Find a step-length lambda using an affine covariant algorithm due to
        % Deuflhard [!!!reference], also described in [!!!thesis]
        
        % Monitors whether we want to accept the current steplength
        accept = 0;
        
        % Indicates whether we are in prediction or correction mode (i. e.
        % whether we are finding the first value of lambda at a given Newton
        % step, or whether we are correcting the value initially predicted for
        % that step).
        initPrediction = 1;
        while ( ~accept )       % Iterate until lambda is accepted
            % Check whether we want to predict a value for lambda. Can only do
            % so once we have taken one Newton step, as it is based on
            % information obtained from the previous step
            if newtonCounter > 0 && initPrediction
                % Compute a prediction value
                mu = (normDeltaOld*nrmDeltaBar)/(mynorm(deltaBar-delta,'fro')*normDelta)*lambda;
                lambda = min(1,mu);
                % Indicate that we will now be in correction mode until next
                % Newton step.
                initPrediction = 0;
            end
            
            if lambda < lambdaMin
                disp('Convergence failure')
                % Take full Newton step
                u = u + delta;
                accept = 1;
                normDeltaOld = normDelta;
                initPrediction = 1;
                lambda = 1;
                continue
            end
            
            % Take a trial step
            uTrial = u + lambda*delta;
            
            uTrialb = uTrial.blocks;
            
            deResFunTrial = N.op(x, uTrialb{:}) - rhs;
            
            % Compute a simplified Newton step, using the current derivative of
            % the operator, but with a new right-hand side.
            [deltaBar, disc] = linsolve(L, deResFunTrial, disc);
            
            % We had two output arguments above, need to negate deltaBar
%             deltaBar = -deltaBar;
            deltaBar = -(L\deResFunTrial);    % Old fashion, to be removed
            % TODO: We also need to update the values of the RHS for the BCs
            % here!
            

            
            % The norm of the simplified Newton step is used to compute a
            % contraction factor
            nrmDeltaBar = mynorm(deltaBar);
            
            % Contraction factor
            cFactor = nrmDeltaBar/normDelta;
            
            muPrime = (.5*normDelta*lambda^2)/...
                (mynorm(deltaBar-(1-lambda)*delta,'fro'));
            
            if cFactor >=1
                lambda = min(muPrime,.5*lambda);
                continue;
            end
            
            lambdaPrime = min(1,muPrime);
            
            if lambdaPrime == 1 && nrmDeltaBar < errTol
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
        
        errEst = NaN;
    else    % We are in undamped phase
        % Update lambda so that we will print correct information
        lambda = 1;
        
        % Take a full Newton step
        u = u + delta;
        
        % Compute a contraction factor and an error estimate. Can only do so
        % once we have taken one step.
        if ( newtonCounter == 0 )
            cFactor = NaN;
        else
            cFactor = normDelta/normDeltaOld;
            
            if cFactor >= 1
                damped = 1; % Have to resort back to damped
                continue
            end
%             contraFactor(newtonCounter-1) = 2*cFactor;
            errEst =  normDelta/(1-cFactor^2);
            errEstDE = errEst;
            errEstVec(newtonCounter) =  errEst;
%             terminate = eEstSob < deltol;
%             omega = 2*nrmDelta/nrmDeltaVec(newtonCounter)^2;
            % TODO: Error estimate
        end
    end
    
    % Evaluate the residual
    ub = u.blocks;
    % Update counter of Newton steps taken
    newtonCounter = newtonCounter +1;
        
    % Store information about the norm of the updates
    normDeltaVec(newtonCounter) = normDelta;
    % Need to store the norm of the current update to use in damping strategy
    normDeltaOld = normDelta;
    
    % Print info to command window, and/or show plot of progress
    displayInfo('iter', u, delta, newtonCounter, normDelta, cFactor, ...
        length(delta{1}), lambda, length(ub{1}), displayFig, displayTimer, pref)
    
    % TODO: Replace with error estimate -- introduce errorTol in cheboppref
%     errEst = normDelta; % .5*(errEstDE + errEstBC)
    if ( errEst < errTol )
        terminate = 1;
    elseif (newtonCounter > maxIter)
        warning('CHEBOP:solvebvpNonlinear','Newton iteration failed.')
        break
    else
        % Linearize around current solution:
        [L, res] = linearize(N, x, ub, []); % flag to negate contraint RHSs.
        
        % Need to subtract RHS from the residual
        res = res - rhs;
        
        % Assign the preferences to the linop.
        L.prefs = pref;
    end
end

% Show final information
displayInfo('final', u, delta, newtonCounter, errEstDE, errEstBC, displayFig, ...
    displayTimer, pref)

% Return useful information in the INFO structure
info.normDelta = normDeltaVec(1:newtonCounter);
info.error = errEst;
end

function out = mynorm(f,type)
if ( isa(f, 'chebmatrix') )
    out = max(cellfun(@(u) get(u, 'vscale'), f.blocks));
else
    out = get(f, 'vscale');
end
end
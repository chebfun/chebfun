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

% TODO: Document properly.

% Store preferences used in the Newton iteration in separate variables
maxIter  = pref.maxIter;
errTol   = pref.errTol;

% Did the user request damped or undamped Newton iteration? Start in mode
% requested (later on, the code can switch between modes).
prefDamped = pref.damped;
damped = prefDamped;

% Minimum allowed steplength
lambdaMin = pref.lambdaMin;

% Assign initial guess to u
u = u0;

% Initialise the independent variable:
x = chebfun(@(x) x, N.domain);

% Print info to command window, and/or show plot of progress
[displayFig, displayTimer] = displayInfo('init', u0, pref);

% Counter for number of Newton steps taken.
newtonCounter = 0;

% Variable that controls whether we want to stop the Newton iteration, either
% because we converged, we have reached the maximum number of iterations, or the
% process has been identified to be nonconvergent. Also initalise related
% control variables.
terminate = 0;
success = 0;
giveUp = 0;
maxIterExceeded = 0;

% Store a vector with information about the norm of the Newton updates
normDeltaVec = zeros(maxIter, 1);

% Initial damping parameter
lambda = 1;

% Need to subtract the rhs from the residual passed in
res = res - rhs;

% Some initializations
dampingInfo.errTol = errTol;
dampingInfo.normDeltaOld = [];
dampingInfo.normDeltaBar = [];
dampingInfo.lambda = lambda;
dampingInfo.lambdaMin = lambdaMin;
dampingInfo.newtonCounter = newtonCounter;
dampingInfo.deltaBar = [];
dampingInfo.damped = damped;
dampingInfo.x = x;

% Start the Newton iteration!
while ( ~terminate )
    
    % Compute a Newton update
    [delta, disc] = linsolve(L, res, pref);
    
    % We had two output arguments above, need to negate delta.
    delta = -delta;

    % Store the norm of the update
    normDelta = N.norm(delta);
    
    dampingInfo.normDelta = normDelta;
    
    % Are we in damped mode?
    if ( damped )
        
        % Find the next Newton iterate (the method finds the step-size, then
        % takes the damped Newton and returns the next iterate)
        [u, dampingInfo] = dampingErrorBased(N, u, rhs, delta, ...
            L, disc, dampingInfo);
        
        % If we're in damped mode, we don't get an error estimate
        errEst = NaN;
        
        % Extract info from the dampingInfo struct, used for printing
        lambda = dampingInfo.lambda;
        cFactor = dampingInfo.cFactor;
        
        % Do we want to keep Newton in damped mode?
        damped = dampingInfo.damped;
        
        % Is the damping strategy telling us to give up?
        giveUp = dampingInfo.giveUp;
        
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
                % Have to resort back to damped (but only if the user wanted
                % damped Newton in the first place).
                damped = prefDamped;
                continue
            end
            
            % Error estimate a la Deuflhard
            errEst =  normDelta/(1-cFactor^2);
        end
    end
    
    % Update counter of Newton steps taken
    newtonCounter = newtonCounter + 1;
    dampingInfo.newtonCounter = newtonCounter;
        
    % Store information about the norm of the updates
    normDeltaVec(newtonCounter) = normDelta;
    % Need to store the norm of the current update to use in damping strategy
    normDeltaOld = normDelta;
    dampingInfo.normDeltaOld = normDeltaOld;
    
    % Grab the blocks of U so that we can print info and evaluate the residual.
    ub = u.blocks;
    
    % Print info to command window, and/or show plot of progress
    displayInfo('iter', u, delta, newtonCounter, normDelta, cFactor, ...
        length(delta{1}), lambda, length(ub{1}), displayFig, displayTimer, pref)
    
    if ( errEst < errTol )  % Sweet, we have converged!      
        success = 1;
    elseif (newtonCounter > maxIter)
        maxIterExceeded = 1;
    else
        % Linearize around current solution:
        [L, res] = linearize(N, ub, x);
        
        % TODO: Ensure this actually takes care of resetting derivative
        % information as well!
        
        % Need to subtract RHS from the residual
        res = res - rhs;
        
        % Assign the preferences to the linop.
        L.prefs = pref;
    end
    
    % Should we stop the Newton iteration?
    if ( success || maxIterExceeded || giveUp )
        terminate = 1;
    end
end

% Evaluate how far off we are from satisfying the boundary conditions.
errEstBC = evalBCnorm(N, u, x);

% Print information depending on why we stopped the Newton iteration
if ( success )

% Show final information.
displayInfo('final', u, delta, newtonCounter, errEst, errEstBC, displayFig, ...
    displayTimer, pref)
elseif ( maxIterExceeded )
    warning('CHEBOP:solvebvpNonlinear:maxIter',...
        ['Newton iteration failed. Maximum number of iterations exceeded.\n',...
        'See help cheboppref for how to increase the number of steps allowed'])
else
    warning('CHEBOP:solvebvpNonlinear:notConvergent',...
        ['Newton iteration failed. Newton iteration is not convergent.\n', ...
        'Please try supplying a better initial guess via the .init field \n' ...
        'of the chebop'])
end

% Return useful information in the INFO structure
info.normDelta = normDeltaVec(1:newtonCounter);
info.error = errEst;
end

function bcNorm = evalBCnorm(N, u, x)
% TODO: This might be useful elsewehere (i.e. chebop/linearize), do we want to
% move this into a separate file?

bcNorm = 0;

uBlocks = u.blocks;

% Evaluate left boundary condition(s):
if ~( isempty(N.lbc) )
    % Evaluate.
    lbcU = N.lbc(uBlocks{:});
    
    % The output might be a CHEBFUN, or a chebmatrix
    if ( isa(lbcU, 'chebfun') )
        bcNorm = bcNorm + sum(feval(lbcU, N.domain(1)).^2);
    elseif ( isa(lbcU, 'chebmatrix') ) 
        % Loop through the components of LBCU.
        for k = 1:numel(lbcU)
            % Obtain the kth element of the CHEBMATRIX
            lbcUk = lbcU{k};
            % Evaluate the function at the left endpoint
            bcNorm = bcNorm + feval(lbcUk, N.domain(1))^2;
        end
    end
end


% Evaluate right boundary condition(s):
if ~( isempty(N.rbc) )
    % Evaluate.
    rbcU = N.rbc(uBlocks{:});
    
    % The output might be a CHEBFUN, or a chebmatrix
    if ( isa(rbcU, 'chebfun') )
        bcNorm = bcNorm + sum(feval(rbcU, N.domain(end)).^2);
    elseif ( isa(rbcU, 'chebmatrix') ) 
        % Loop through the components of RBCU.
        for k = 1:numel(rbcU)
            % Obtain the kth element of the CHEBMATRIX
            rbcUk = rbcU{k};
            % Evaluate the function at the left endpoint
            bcNorm = bcNorm + feval(rbcUk, N.domain(1))^2;
        end
    end
end


% Evaluate and linearise the remaining constraints:
if ( ~isempty(N.bc) )
    % Evaluate. The output, BCU, will be a vector.
    bcU = N.bc(x, uBlocks{:});
    
    bcNorm = bcNorm + sum(norm(bcU).^2);
end

bcNorm = sqrt(bcNorm);

end

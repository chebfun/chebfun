function [u, info] = solvebvpNonlinear(N, rhs, L, u0, res, pref, displayInfo)
%SOLVEBVPNONLINEAR      Solve a nonlinear BVP, using damped Newton iteration.
% The inputs to the method are:
%   N:      Nonlinear CHEBOP
%   rhs:    A CHEBMATRIX, right hand side of ODE
%   L:      A LINOP, linearisation of N around the initial guess
%   u0:     A CHEBMATRIX, an initial guess of solution
%   res:    A CHEBMATRIX, residual of ODE at initial guess
%   pref:   CHEBOPPREF preference structure
%
% The outputs are
%   u:      A CHEBMATRIX, that represents the solution of the BVP if the method
%           converged to a solution.
%   info:   A MATLAB struct with useful information about the Newton iteration.
%
% See also: CHEBOP/SOLVEBVP, CHEBOP/DAMPINGERRORBASED.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   This method also accepts the function handle DISPLAYINFO, that allows
%   passing in a function handle to a displaying method that is called during
%   the Newton iteration. This allows separating the displaying process for
%   regular CHEBOP use and CHEBGUI. See chebop/displayInfo() and
%   chebgui/displayInfo() for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store preferences used in the Newton iteration in separate variables
maxIter  = pref.maxIter;
errTol   = pref.errTol;

% Did the user request damped or undamped Newton iteration? Start in mode
% requested (later on, the code can switch between modes).
prefDamping = pref.damping;
damping =      prefDamping;

% Assign initial guess to u:
u = u0;

% Initialise the independent variable:
x = chebfun(@(x) x, N.domain);

% Counter for number of Newton steps taken.
newtonCounter = 0;

% Variables that controls whether we want to stop the Newton iteration, either
% because we converged, we have reached the maximum number of iterations, or the
% process has been identified to be nonconvergent.
success = 0;
giveUp = 0;
maxIterExceeded = 0;
terminate = 0;

% Store a vector with information about the norm of the Newton updates:
normDeltaVec = zeros(maxIter, 1);

% Initial damping parameter:
lambda = 1;

% Need to subtract the rhs from the residual passed in.
res = res - rhs;

% Initial estimate of error
errEst = inf;

% Some initializations of the DAMPINGINFO struct. See 
%   >> help dampingErrorBased 
% for discussion of this struct. 
dampingInfo.errTol =        errTol;
dampingInfo.normDeltaOld =  [];
dampingInfo.normDeltaBar =  [];
dampingInfo.lambda =        lambda;
dampingInfo.lambdaMin =     pref.lambdaMin;
dampingInfo.newtonCounter = newtonCounter;
dampingInfo.deltaBar =      [];
dampingInfo.damping =       damping;
dampingInfo.x =             x;
dampingInfo.giveUp =        0;

linpref = pref;
linpref.errTol = pref.errTol/10;

% Get the differential order of the LINOP L (needed when evaluating the residual
% of periodic boundary conditions):
diffOrder = L.diffOrder;

% Start the Newton iteration!
while ( ~terminate )
    
    % Compute a Newton update:
    [delta, disc] = linsolve(L, res, linpref, vscale(u));

    % We had two output arguments above, need to negate DELTA.
    delta = -delta;

    % Store the norm of the update:
    normDelta = norm(delta);
    
    % Assign to the DAMPINGINFO struct:
    dampingInfo.normDelta = normDelta;
    
    % At the first Newton iteration, we have to do additional checks.
    if ( newtonCounter == 0)
        % Did we actually get an initial passed that solves the BVP?
        if ( normDelta/sum(vscale(u)) < errTol/100 )
            displayInfo('exactInitial', pref);
            info.error = NaN;
            info.normDelta = normDelta;
            return
        else
            % We actually have to start the Newton iteration. Print info to
            % command window, and/or show plot of progress:
            [displayFig, displayTimer] = displayInfo('init', u0, pref);
        end
    end
    
    % Are we in damped mode?
    if ( damping )
        
        % Find the next Newton iterate (the method finds the step-size, then
        % takes the damped Newton and returns the next iterate).
        [u, dampingInfo] = dampingErrorBased(N, u, rhs, delta, ...
            L, disc, dampingInfo);
        
        % If we're in damped mode, we don't get an error estimate...
        errEst = NaN;
        
        % Extract info from the dampingInfo struct, used for printing:
        lambda = dampingInfo.lambda;
        cFactor = dampingInfo.cFactor;
        
        % Do we want to keep Newton in damped mode?
        damping = dampingInfo.damping;
        
        % Is the damping strategy telling us to give up?
        giveUp = dampingInfo.giveUp;
        
        % Did we converge within the damped phase?
        success = dampingInfo.success;
        
    else % We are in undamped phase
        
        % Update lambda so that we will print correct information in the
        % displayInfo() method.
        lambda = 1;
        
        % Take a full Newton step:
        u = u + delta;
        
        % Compute a contraction factor and an error estimate. Can only do so
        % once we have taken one step.
        if ( newtonCounter == 0 )
            cFactor = NaN;
        else
            % Compute the contraction factor of this iterate:
            cFactor = normDelta / normDeltaOld;
            
            if ( cFactor >= 1 )
                % We're not observing the nice convergence of Newton iteration
                % anymore. Have to resort back to damped iteration (but only if
                % the user wanted damped Newton in the first place).
                damping = prefDamping;
                if ( damping ) 
                    continue    % Go back to the start of loop
                end
            else
                % Error estimate based on the norm of the update and the contraction
                % factor.
                errEst =  normDelta / (1 - cFactor^2);
            end
        end
        
    end
    
    % Update counter of Newton steps taken:
    newtonCounter = newtonCounter + 1;
    dampingInfo.newtonCounter = newtonCounter;
        
    % Store information about the norm of the updates:
    normDeltaVec(newtonCounter) = normDelta;
    
    % Need to store the norm of the current update to use in damping strategy:
    normDeltaOld = normDelta;
    dampingInfo.normDeltaOld = normDeltaOld;
    
    % Find the maximum length of the current solution:
    len = max(cellfun(@length, u.blocks(:)));
    
    % Print info to command window, and/or show plot of progress
    [displayTimer, stopReq] = displayInfo('iter', u, delta, newtonCounter, ...
        normDelta, cFactor, length(delta{1}), lambda, len, displayFig, ...
        displayTimer, pref);
    
    if ( errEst < errTol )  
        % Sweet, we have converged!      
        success = 1;
    elseif ( newtonCounter > maxIter )
        % Damn, we failed.
        maxIterExceeded = 1;
    elseif ( stopReq )
        % User requested to stop the iteration. Generally, this will only happen
        % in GUI mode.
        
        % Do nothing.
    else
        % Linearize around current solution:
        [L, res] = linearize(N, u, x);
        % Need to subtract the original RHS from the residual:
        res = res - rhs;
    end
    
    % Should we stop the Newton iteration?
    if ( success || maxIterExceeded || ( giveUp == 1) || stopReq )
        break
    end
    
end

% Simplify the result before returning it and printing solver info:
u = simplify(u);

% Evaluate how far off we are from satisfying the boundary conditions.
errEstBC = normBCres(N, u, x, diffOrder, pref);

% Print information depending on why we stopped the Newton iteration.
if ( success )
    % Show final information.
    displayInfo('final', u, delta, newtonCounter, errEst, errEstBC, ...
        displayFig, displayTimer, pref)
elseif ( maxIterExceeded )
    warning('CHEBFUN:CHEBOP:solvebvpNonlinear:maxIter',...
        ['Newton iteration failed. Maximum number of iterations exceeded.\n',...
        'See help cheboppref for how to increase the number of steps allowed.'])
elseif ( giveUp )
    warning('CHEBFUN:CHEBOP:solvebvpNonlinear:notConvergent',...
        ['Newton iteration failed.\n', ...
        'Please try supplying a better initial guess via the .init field \n' ...
        'of the chebop.'])
end

% Return useful information in the INFO structure
info.normDelta = normDeltaVec(1:newtonCounter);
info.error = errEst;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bcNorm = normBCres(N, u, x, diffOrder, pref)
%NORMBCRES   Compute residual norm of the boundary conditions.
%   NORMBCRES(N, U, X, DIFFORDER, PREF) returns the combined Frobenius norm of N.lbc(U),
%   N.rbc(U), and N.bc(X, U).

% [TODO]: This might be useful elsewehere (i.e. chebop/linearize), do we want to
% move this into a separate file?

% Initialize:
bcNorm = 0;

% If nargin(N) <= 2, but the dimension of the solution guess passed is greater
% than 1, we are working with the @(x,u) [diff(u{1}) + u{2}; ...] syntax. Need
% to make the code aware of this.
if ( nargin(N) <= 2 && max(size(u)) > 1 )
    cellArg = 1;
else
    cellArg = 0;
end

% Extract the blocks from the CHEBMATRIX U.
uBlocks = u.blocks;

% Evaluate left boundary condition(s):
if ( ~isempty(N.lbc) )
    % Evaluate.
    if ( cellArg )
        lbcU = N.lbc(u);
    else
        lbcU = N.lbc(uBlocks{:});
    end
    
    % The output might be a CHEBFUN, or a CHEBMATRIX
    if ( isa(lbcU, 'chebfun') )
        bcNorm = bcNorm + sum(feval(lbcU, N.domain(1)).^2);
    elseif ( isa(lbcU, 'chebmatrix') ) 
        % Loop through the components of LBCU.
        for k = 1:max(size(lbcU))
            % Obtain the kth element of the CHEBMATRIX
            lbcUk = lbcU{k};
            % Evaluate the function at the left endpoint
            bcNorm = bcNorm + feval(lbcUk, N.domain(1))^2;
        end
    end
end

% Evaluate right boundary condition(s):
if ( ~isempty(N.rbc) )
    % Evaluate.
    if ( cellArg )
        rbcU = N.rbc(u);
    else
        rbcU = N.rbc(uBlocks{:});
    end
    
    % The output might be a CHEBFUN, or a CHEBMATRIX
    if ( isa(rbcU, 'chebfun') )
        bcNorm = bcNorm + sum(feval(rbcU, N.domain(end)).^2);
    elseif ( isa(rbcU, 'chebmatrix') ) 
        % Loop through the components of RBCU.
        for k = 1:max(size(rbcU))
            % Obtain the kth element of the CHEBMATRIX
            rbcUk = rbcU{k};
            % Evaluate the function at the left endpoint
            bcNorm = bcNorm + feval(rbcUk, N.domain(end))^2;
        end
    end
end

% Evaluate and linearise the remaining constraints:
disc = pref.discretization();
tech = disc.returnTech();
techUsed = tech();

if ( ~isempty(N.bc) || isequal(pref.discretization, @trigcolloc) )
    % Periodic case. 
    if ( (isa(N.bc, 'char') && strcmpi(N.bc, 'periodic')) || ...
            isPeriodicTech(techUsed) )
        bcU = 0;
        % Need to evaluate the residual of the boundary condition for each
        % independent variable uBlocks{k} separately, since each variable can
        % have a different maximum differential order associated with it in a
        % problem.
        for k = 1:numel(uBlocks)
            for l = 0:max(diffOrder(:, k))
                % Compute residual of appropriately many derivatives:
                bcU = bcU + (feval(diff(uBlocks{k}, l), N.domain(end)) - ...
                    feval(diff(uBlocks{k}, l), N.domain(1)))^2;
            end
        end
        bcNorm = bcNorm + bcU;
        
    else
        % Evaluate. The output, BCU, will be a vector.
        bcU = N.bc(x, uBlocks{:});
        bcNorm = bcNorm + norm(bcU, 2).^2;
    end
end

bcNorm = sqrt(bcNorm);

end

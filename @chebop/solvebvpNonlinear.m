function u = solvebvpNonlinear(N, rhs, pref)

% Store preferences used in the Newton iteration in separate variables
maxIter  = pref.maxIter;
discType = pref.discretization;
delTol   = pref.deltol;
% plotMode determines whether we want to stop between plotting iterations
plotMode = lower(pref.plotting);
% Did the user request damped or undamped Newton iteration?
damped = pref.damped;

% NUMVARS indicate how many unknown function we seek.
numVars = nargin(N.op) - 1;

% Store the domain we're working with.
dom = N.domain;

% Initialise a zero ADCHEBFUN:
zeroFun = chebfun(0, dom);
u0 = cell(numVars, 1);
for k = 1:numVars
    u0{k} = zeroFun;
end
u0 = chebmatrix(u0);

% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Print info to command window, and/or show plot of progress
[displayFig, displayTimer] = N.displayInfoInit(u0, pref);

% Linearise
[L, affine, isLinear] = linearise(N, x, u0);

    function out = mynorm(f)
        if ( isa(f, 'chebmatrix') )
            out = max(cellfun(@(u) get(u, 'vscale'), f.blocks));
        else
            out = get(f, 'vscale');
        end
    end

if ( ~isempty(N.init) )
    u = N.init;
    L = linearise(N, x, u);
    res = N.op(x, u{:}) - rhs;
    L.discretization = discType;
    du = L\res;
    %                 else
    %                     u = makeGuess(L, x);
    %                     L = linearise(N, x, u);
    %                     res = N.op(x, u{:}) - rhs;
    %                     du = L\res;
else
    u = u0;
    L.discretization = discType;
    du = L\(rhs - affine);
end

u = u - du;
ub = u.blocks;
res = N.op(x, ub{:}) - rhs;

% fprintf('step\t normUpdate\t\t  normRes\n')
normUpdate(1,1) = mynorm(du);
normRes(1,1) = mynorm(res);
% fprintf(' %2.2d\t%16.16f\t %16.16f\n', 1, normUpdate(1,1), normRes(1,1))

% Counter for number of Newton steps taken.
newtonCounter = 0;

% Variable that controls whether we want to stop the Newton iteration, either
% because we converged, or the process has been identified to be nonconvergent.
terminate = 0;
while ( ~terminate ) && ( normRes(end) > 1e-12)
    % Linearise around current solution:
    L = linearise(N, x, ub, []); % flag to negate contraint RHSs.
    % Solve the linearised system:
    L.discretization = discType;
    du = L\res;
    % Append the Newton step:
    u = u - du;
    % Evaluate the residual
    ub = u.blocks;
    res = N.op(x, ub{:}) - rhs;
    
    % Update counter of Newton steps taken
    newtonCounter = newtonCounter +1;
    
    % Stop if well converged, or stagnated:
    normDelta = mynorm(du);
    normUpdate(newtonCounter,1) = normDelta;
    normRes(newtonCounter,1) = mynorm(res);
    
    % Print info to command window, and/or show plot of progress
    N.displayInfoIter(u, newtonCounter, normDelta, length(du{1}), ...
        displayFig, displayTimer, pref)

    if ( normUpdate(newtonCounter,1) < 1e-12 )
        terminate = 1;
        %                     elseif (newt > 3 && normUpdate(newt) > 0.1*normUpdate(newt-3))
        %                         warning('CHEBFUN:bvpsc','Newton iteration stagnated.')
        %                         break
    elseif (newtonCounter > maxIter)
        warning('CHEBFUN:bvpsc','Newton iteration failed.')
        break
    end
    

    
end

end
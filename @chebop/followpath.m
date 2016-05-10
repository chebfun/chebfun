function [uquasi, lamvec, mvec, lamfun, mfun] = followpath(N, lam0, varargin)
%FOLLOWPATH    A pseudo-arclength continuation algorithm for ODEs in Chebfun.
%
% Calling sequence:
%   [U, LAM] = FOLLOWPATH(N, LAM0, 'OPT1', VAL1, ...)
%
% Here, the inputs are:
%   
%   N    : A chebop, whose N.op arguments are x, u and lambda, and boundary
%          conditions also depend on u and lambda.
%   LAM0 : Initial value of lambda for finding an initial point on the solution
%          curve.
%
% It is possible to pass the method various option pairs on the form
%   'OPTIONNAME', OPTIONVALUE
% The options supported are
%   'UINIT'     : Initial solution U on the solution curve, with LAMBDA = LAM0.
%   'DIRECTION' : Whether the curve should be tracked in positive or negative
%                 direction. Possible values: +1 (default), -1.
%   'MEASURE'   : An anonymous function that takes U as argument, used for
%                 computing data for a bifurcation diagram during the tracing of
%                 the solution curve. See note below on how to call the method
%                 for the values of MEASURE to be returned.
%   'PLOTTING'  : Whether intermediate solutions and bifurcation diagram should
%                 be plotted during the tracing of the curve. Possible values
%                 TRUE, FALSE (default). Note that if 'PLOTTING' is set to TRUE,
%                 'MEASURE' has to be passed as well.
%   'MAXSTEPNO' : Maximum number of points to compute on the solution curve.
%                 Default value: 25.
%   'STEPMAX'   : Maximum steplength accepted for pseudo-arclength continuation.
%                 Default value: .5.
%   'STEPMIN'   : Minimum steplength accepted for pseudo-arclength continuation.
%                 Note that if the path-following algorithm wants to take a step
%                 smaller than STEPMIN, the program will terminate.
%                 Default value: 1e-4.
%   'STEPINIT'  : Initial steplength accepted for pseudo-arclength continuation.
%                 Default value: Set to be the same as STEPMAX.
%   'STOPFUN'   : An anonymous function that takes U and LAMBDA as arguments and
%                 returns a Boolean value, so that when 
%                     STOPFUN(U,LAMBDA) == TRUE
%                 the pathfollowing program gets terminated, even if MAXSTEPNO
%                 has not been reached.
%   'PREFS'     : A CHEBOPPREF object. By default, a new CHEBOPPREF object with
%                 the current global preferences is used if PREFS is not passed.
%                 Note that the discretization option of PREFS has to be
%                 specified as a function handle, not string, see 'help
%                 cheboppref' for more details.
%
% The outputs are
%   U   : An array-valued CHEBFUN that contains all computed functions on the
%         solution curve.
%   LAM : A vector containing all values of the parameter LAMBA computed on the
%         solution curve, so that the ith element of LAM corresponds to the ith
%         column of U.
%
% Note 1: If no UINIT is passed, the initial solution U used is the one computed
% by the CHEBOP SOLVEBVP algorithm.
%
% Note 2: If MEASURE is passed, it is possible to call the method with three
% outputs:
%   [UQUASI, LAMVEC, MVEC] = FOLLOWPATH(N, LAM0, 'OPT1', VAL1, ...)
% In this case, MVEC contains the value of MEASURE at all points on the solution
% curve computed.
%
% Note 3: Calling the function with five outputs:
%   [UQUASI, LAMVEC, MVEC, LAMFUN, MFUN] = FOLLOWPATH(N, LAM0, ...)
% also returns the chebfuns LAMFUN and MFUN, which are spline interpolants of
% the LAMVEC and MVEC data. For a smooth bifurcation diagram, it is then
% possible to call
%   plot(LAMFUN, MFUN)
%
% Example 1 -- Bratu problem, continuation on lambda parameter:
%   N = chebop(@(x,u,lam) diff(u,2) + lam*exp(u), [0 1]);
%   N.lbc = @(u,lam) u;
%   N.rbc = @(u,lam) u;
%   lam0 = 0.01;
%   % Call method, no plotting, no printing
%   [u, lamvec] = followpath(N, lam0);
%   % Call method, specifying more options
%   [u, lamvec, mvec] = followpath(N, lam0, ...
%       'measure', @(u) u(.5), 'printing', true, 'plotting',true);
%   
% Example 2 -- Herceg problem (singularly perturbed ODE). Fix solution value at
%              left endpoint, vary slope at right endpoint.
%   d = [0 1];
%   ep = 2^-2;
%   % Start by finding initial solution U on curve
%   N = chebop(@(x,u) -ep^2*diff(u,2)+(u.^2+u-.75).*(u.^2+u-3.75), d);
%   N.lbc = 0;
%   N.rbc = 0;
%   u0 = N\0;
%   % Trace solution curve:
%   N = chebop(@(x,u,lam) -ep^2*diff(u,2) + (u.^2+u-.75).*(u.^2+u-3.75), d);
%   lam0 = feval(diff(u0),d(2)); % Initial value for LAMBDA
%   N.lbc = @(u, lam) u;
%   N.rbc = @(u, lam) diff(u) - lam;
%   [u, lamvec, mvec, lamfun, mfun] = followpath(N, lam0, 'maxstepno', 30, ...
%       'uinit', u0, 'measure', @(u)u(1), 'stepmax', 1, 'printing', 1);
%   % Plot a bifurcation diagram
%   figure, plot(lamfun, mfun)
%
% Example 3 -- Herceg problem, continuation on perturbation parameter
%   ep = 2^-2;
%   H = chebop(@(x,u,lam) -lam*diff(u,2)+(u.^2+u-.75).*(u.^2+u-3.75), [0 1]);
%   lam0 = ep^2;
%   H.lbc = @(u, lam) u;
%   H.rbc = @(u, lam) u;
%   measure = @(u) norm(diff(u), inf); 
%   [u, lamvec, mvec] = followpath(H, lam0, ...
%       'measure', measure, 'direction', -1, 'plotting', 1, ...
%       'stopfun', @(u,lam) lam < 5e-3, 'stepmax', .1);

% Copyright 2016 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% ========================= Developer comment ==================================
% For every BVP we solve during the run of FOLLOWPATH, we need to compute
% derivatives of the operators involved. Currently this is done with AD, which
% is fairly slow as it has to happen so often. It'd be desirable to be able to
% pass in a recipe for the derivative as an option, e.g.:
%   'DER'       : A struct that contains anonymous functions that describe the
%                 Frechet derivatives of N. Passing DER leads to a faster run of
%                 the program, as then it is not necessary to user automatic
%                 differentiation to compute the operator derivatives required.
% or analyse the evaluation tree of the operator to compute the derivative once
% and for all at the start of the run.

% Set default values
uinit = [];                 % By default, use chebop to compute the initial sol
measure = [];               % By default, no measure is computed during run
plotting = false;           % By default, don't plot anything during the run
printing = false;           % By default, don't print anything during the run
direction = 1;              % By default, go in positive direction
stepmax = .5;               % Maximum steplength
stepmin = 1e-4;             % Mininum steplength
stepinit = [];              % Initial steplength
maxstepno = 25;             % Maximum number of steps taken
stopfun = @(u, lambda) 0;   % Default stopfun that always evaluates to false
prefs = [];


% Parse VARARGIN
while ( ~isempty(varargin) ) % Go through all elements
    if ( ~ischar(varargin{1}) && ~isnumeric(varargin{2}) )
        error('followpath:inputArgument','Incorrect options input arguments');
    end
    val = varargin{2};
    switch lower(varargin{1})
        case 'uinit'
            uinit = val;
        case 'direction'
            direction = val;
        case 'measure'
            measure = val;
        case 'plotting'
            plotting = val;
        case 'maxstepno'
            maxstepno = val;
        case 'stepmax'
            stepmax = val;
        case 'stepmin'
            stepmin = val;
        case 'stepinit'
            stepinit = val;
        case 'printing'
            printing = val;
        case 'stopfun'
            stopfun = val;
        case 'prefs'
            prefs = val;
    end
    
    % Throw away option name and argument and move on
    varargin(1:2) = [];
end

% If STEPINIT was not specified, set it to STEPMAX
if ( isempty(stepinit) )
    stepinit = stepmax;
end


% No initial u passed => find a u0 matching lam0:
if ( isempty(uinit) )
    
    if ( printing )
        fprintf('=====================================================\n')
        fprintf('Computing initial solution for pathfollowing...')
    end
    
    % Create an "initial chebop" from the augmented one passed in. 
    Ninit = N;
    Ninit.op = @(x,u) N.op(x, u, lam0);
    % This requires us to evaluate the boundary conditions, so that we get BCs
    % that only depend on U (by fixing LAM to be LAM0):
    if ( ~isempty(Ninit.bc) )
        Ninit.bc = @(x,u) Ninit.bc(x,u,lam0);
    end
    if ( ~isempty(Ninit.lbc) )
        Ninit.lbc = @(u) Ninit.lbc(u,lam0);
    end
    if ( ~isempty(Ninit.rbc) )
        Ninit.rbc = @(u) Ninit.rbc(u,lam0);
    end
    
    % Compute the initial solution
    uinit = Ninit\0;
    
    if ( printing )
        fprintf(' done!\n')
        fprintf('=====================================================\n\n')
    end
end

% No preferences specified, use current CHEBOPPREF
if ( isempty(prefs) )
    prefs = cheboppref();
    
    % By default, cheboppref now has the discretization option specified as a
    % string ('values' or 'coeffs'). As we're not going through CHEBOP/SOLVEBVP
    % below, we need to call determineDiscretization to get the correct
    % discretization option as a function handle. That method require a linop
    % argument as well (in case of breakpoints), hence the call to linearize:
    L = linearize(N, [uinit; lam0]);
    lengthDom = max(length(L.domain), length(uinit.domain));
    prefs = determineDiscretization(N, lengthDom, prefs);
end

% Did we have a measure passed?
if ( isempty(measure) )
    haveMeasure = false;
else
    haveMeasure = true;
end

% Store all the solutions to be returned, as well as lambda values and the
% measure values. The sizes are all maxstepno + 1, as we also store the initial
% solution.
uquasi = cell(1, maxstepno + 1);
uquasi{1} = uinit;
lamvec = zeros(maxstepno + 1, 1);
lamvec(1) = lam0;
mvec = zeros(maxstepno + 1, 1);

% If we have a measure passed, we evaluate it at the initial solution, and
% convert it to a nice string for the plot:
if ( haveMeasure )
    measu = measure(uinit);
    mvec(1) = measu;
    % Obtain a nice string to set on ylabel of bifurcation diagram
    mstring = func2str(measure);
    % Throw away the @(u) part
    mstring = mstring(min(strfind(mstring,')'))+1:end);
end

% If user wants to plot, a measure has to be passed!
if ( plotting && ~haveMeasure )
    error('CHEBFUN:CHEBOP:followpath:plotButNoMeasure', ...
        'If plotting is ON for path-following, measure has to be supplied.')
end

% If plotting, before starting path following, plot initial information
if ( plotting )
    % Plot solution
    subplot(1,2,1);
    plot(uinit)
    title(['Solution for \lambda =' num2str(lam0)])
    xlabel('x'),ylabel('u(x)')
    
    % Plot the bifurcation diagram
    subplot(1,2,2)
    plot(lamvec,mvec,'-*')
    title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
    drawnow, shg
end

% Create the independent problem on the domain of the problem:
dom = domain(uinit);
x = chebfun(@(x) x, domain(uinit));
% A diagonal sum operator on the domain of the problem, used below.
dSum = diag(sum(dom));

% Set up initial tangent and tau we want to be orthogonal to.
told = chebfun(0, domain(uinit));
tauold = 1;

% Print initial information if we're printing:
if ( printing )
    if ( haveMeasure )
        fprintf('#Sol    #Newton     lambda     Steplength    Measure    \n')
        fprintf('--------------------------------------------------------\n')
        fprintf('%3i \t       \t   %6.2e \t            %6.2e \n', ...
            1, lam0, measu)
    else
        fprintf('#Sol    #Newton     lambda     Steplength   \n')
        fprintf('--------------------------------------------\n')
        fprintf('%3i \t       \t   %6.2e \t \n', ...
            1, lam0)
    end    
end

% Variable for keeping track of whether we accept the tangent step or retract.
% If retract == 1, the Newton step told us to go back along the tangent and
% shrink the stepsize.
retract = false;
% Counter for the number of steps we take:
counter = 1;
% Store previous solution:
uold = uinit;
% Store previous value of lambda:
lamold = lam0;
% Set steplength to initial steplength specified:
sl = stepinit;

% Start the pseudo-arclength path following algorithm. It proceeds as follows:
%   1. Find a tangent direction.
%   2. Move in the direction of the tangent for the given steplength.
%   3. Compute the Newton correction.
%   4. If Newton was happy, accept the new point on the curve. If not, shrink
%      the steplength by a factor of 4 and go back to step 2.
while ( counter < maxstepno )
    % Find a tangent direction, but only if we were told by Newton not to
    % retract. At the start of the while loop, recall that RETRACT == false.
    if ( ~retract )
        % Compute the tangent
        [t, tau] = tangentBVP(N, uold, lamold, told, tauold, x, dSum, prefs);
        % At the start, we need to ensure we're going in right direction:
        if ( counter == 1 )
            t = direction*t;
            tau = direction*tau;
        end
        % Move in the direction of the tangent
        uinit = uold + sl*t;
        laminit = lamold + sl*tau;
    end

    % Find a Newton correction to get back on the solution curve:
    [u, lam, newtonIter, retract] = ...
        newtonBVP(N, uinit, laminit, t, tau, x, dSum, prefs);
    
    % If the Newton correction algorithm told us we were trying to take too long
    % tangent steps, we decrease the steplength.
    if ( retract )
        % Move in the direction of the current tangent, but only with
        % quarter of the steplength
        sl = sl/4;
        uinit = uold + sl*t;
        laminit = lamold + sl*tau;
        
        % Have we reached the mininum steplength approved?
        if ( sl < stepmin )
            disp('FAILED: sl < stepmin')
            break
        end
        
        % Go back to the start of the while loop, haven taken a smaller tangent
        % step:
        continue
    end
    
    % We've found a new point, update counter:
    counter = counter + 1;
    
    % Store values for plotting
    if ( haveMeasure )
        measu = measure(u);
        mvec(counter) = measu;
    end
    lamvec(counter) = lam;
    
    % Print information at the current step.
    if ( printing )
        if ( haveMeasure )
            fprintf('%3i \t   %2i \t   %6.2e \t %6.4f     %6.2e \n', ...
                counter, newtonIter, lam, sl, measu)
        else
            fprintf('%3i \t   %2i \t   %6.2e \t %6.4f \n', ...
                counter, newtonIter, lam, sl)
        end
    end

    % Update diagram if we're plotting:
    if ( plotting )
        % Plot current solution
        subplot(1,2,1);
        plot(u)
        title(['Solution for \lambda =' num2str(lam)])
        xlabel('x'),ylabel('u(x)')
        
        % Update bifurcation diagram:
        subplot(1,2,2)
        % Create splines to draw a smooth bifurcation curve
        lamspline = chebfun.spline(linspace(0, 1, counter), ...
            lamvec(1:counter));
        mspline = chebfun.spline(linspace(0, 1, counter), ...
            mvec(1:counter));
        plot(lamspline, mspline)
        
        % Add the points on the solution curve
        hold on
        set(gca, 'ColorOrderIndex', 1)
        plot(lamvec, mvec, '*')
        hold off
        title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
        drawnow, shg
    end
    
    % If we're experiencing good Newton convergence, we try to get the
    % steplength closer to the maximum steplength allowed:
    if ( newtonIter <= 3 )
        sl = min(sl*2, stepmax);
    end
    
    % If we've been successful to get here, update old variable values:
    told = t;
    tauold = tau;
    uold = u;
    lamold = lam;
    
    % Update quasimatrix to be returned
    uquasi{counter} = u;
    
    % Is STOPFUN telling us to stop?
    if ( stopfun(u, lam) )
        break
    end
    
end

% Throw away unneeded elements of the cell and vectors:
uquasi(counter + 1 : end) = [];
lamvec(counter + 1 : end) = [];
mvec(counter + 1 : end) = [];

% Convert cell to an array valued chebfun to be returned:
uquasi = chebfun(chebmatrix(uquasi));

% Create splines of lambdas and measures to be returned for plotting afterwards:
lamfun = chebfun.spline(linspace(0, 1, length(lamvec)),lamvec);
mfun = chebfun.spline(linspace(0, 1, length(lamvec)), mvec);

end

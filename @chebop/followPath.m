function [uquasi, lamvec, mvec] = followPath(N, lam0, varargin)
%FOLLOWPATH    A pseudo-arclength continuation algorithm for ODEs in Chebfun
%
% Calling sequence:
%   [UQUASI, LAMVEC] = FOLLOWPATH(N, LAM0, 'OPT1', VAL1, ...)
% Here, the inputs are:
%   
%   N    : A chebop, whose N.op arguments are x, u and lambda.
%   lam0 : Initial value of lambda for finding initial point on solution curve.
%
% It is possible to pass the method various option pairs on the form
%   'OPTIONNAME', OPTIONVALUE
% The options supported are
%   'UINIT'     : Initial solution U on the solution curve.
%   'DIRECTION' : Whether the curve should be tracked in positive or negative
%                 direction. Possible values: +1 (default), -1.
%   'MEASURE'   : An anonymous function that takes U as argument, used for
%                 drawing a bifurcation diagram during the tracing of the
%                 solution curve.
%   'PLOTTING'  : Whether intermediate solutions and bifurcation diagram should
%                 be plotted during the tracing of the curve. Possible values
%                 TRUE, FALSE (default). Note that if 'PLOTTING' is set to TRUE,
%                 'MEASURE' has to be passed as well.
%   'MAXITER'   : Maximum number of points to compute on the solution curve.
%                 Default value: 25.
%   'SMAX'      : Maximum steplength accepted for pseudo-arclength continuation.
%                 Default value: .5.
%   'SMIN'      : Minimum steplength accepted for pseudo-arclength continuation.
%                 Note that if the path-following algorithm wants to take a step
%                 smaller than SMIN, the program will give up.
%                 Default value: 1e-4.
%   'S0'        : Initial steplength accepted for pseudo-arclength continuation.
%                 Default value: Set to be the same as SMAX.
%   'STOPFUN'   : An anonymous function that takes U and LAMBDA as arguments and
%                 returns a Boolean value, so that when 
%                     STOPFUN(U,LAMBDA) == TRUE
%                 the pathfollowing program gets terminated, even if MAXITER has
%                 not been reached.
%   'DER'       : A struct that contains anonymous functions that describe the
%                 Frechet derivatives of N. [TODO: Describe further and
%                 implement].
%   'PREFS'     : A CHEBOPPREF object. By default, a new CHEBOPPREF object with
%                 the current global preferences is used.
%
% The outputs are
%   U   : An array-valued CHEBFUN that contains all computed functions on the
%         solution curve.
%   LAM : A vector containing all values of the parameter LAMBA computed on the
%         solution curve, so that the ith element of LAM corresponds to the ith
%         column of U.
%
% Note 1: If no UNIT is passed, the initial solution U used is the one computed
% by chebop.
%
% Note 2: If MEASURE is passed, it is possible to call the method with three
% outputs:
%   [UQUASI, LAMVEC, MVEC] = FOLLOWPATH(N, LAM0, 'OPT1', VAL1, ...)
% In this case, MVEC contains the value of MEASURE at all points on the solution
% curve computed.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Set default values
uinit = [];
measure = [];
plotting = false; % Option for plotting
printing = false;
direction = 1;
smax = .5;      % Maximum steplength
smin = 1e-4;    % Mininum steplength
s0 = [];        % Initial steplength
maxiter = 25;
stopfun = @(u, lambda) 0;
prefs = [];

if ( nargin < 8 )
    direction = 1;
end

% Parse varargin
while ~isempty(varargin)  % Recurse
    if ~ischar(varargin{1}) && ~isnumeric(varargin{2})
        error('followpath:inputArgument','Incorrect input arguments');
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
        case 's0'
            s0 = val;
        case 'maxiter'
            maxiter = val;
        case 'slmax'
            smax = val;
        case 'printing'
            printing = val;
        case 'stopfun'
            stopfun = val;
        case 'prefs'
            prefs = val;
    end
    
    % Throw away arguments and move on
    varargin(1:2) = [];
end

% If S0 was not specified, set it to SMAX
if ( isempty(s0) )
    s0 = smax;
end

% No preferences specified
if ( isempty(prefs) )
    prefs = cheboppref();
    % TODO: Deal with values/coeffs, and periodic case.
    prefs.discretization = @chebcolloc2;
end

% No initial u passed -> find a u0 matching lam0
if ( isempty(uinit) )
    if ( printing )
        fprintf('=====================================================\n')
        fprintf('Computing initial solution for pathfollowing...')
    end
    Ninit = N;
    Ninit.op = @(x,u) N.op(x, u, lam0);
    uinit = Ninit\0;
    
    if ( printing )
        fprintf(' done!\n')
        fprintf('=====================================================\n\n')
    end
end

% Did we have a measure passed?
if ( isempty(measure) )
    haveMeasure = false;
else
    haveMeasure = true;
end

% Store all the solutions to be returned
uquasi = uinit;

% Constraint for tangent
J = @(u,lam) sum(u).^2+lam.^2;

% Iterate along path.
% Begin by finding a tangent direction, then set steplength, then compute
% Newton correction, and repeat.
counter = 1;
uold = uinit;
lamold = lam0;
if isa(lam0,'chebconst')
    lamvec = lam0.vals;
else
    lamvec = lam0;
end

if ( haveMeasure )
    mvec = measure(uinit);
    
    % Obtain a nice string to set on ylabel of bifurcation diagram
    mstring = func2str(measure);
    mstring = mstring(min(strfind(mstring,')'))+1:end); % Throw away the @(u) part
end

% If user wants to plot, a measure has to be passed!
if ( plotting && ~haveMeasure )
    error('CHEBFUN:CHEBOP:followPath:plotButNoMeasure', ...
        'If plotting is ON for path-following, measure has to be supplied.')
end

sl = s0;

% If plotting, before starting path following, plot initial information
if ( plotting )
    subplot(1,2,1);
    plot(uinit), title(['Solution for \lambda =' num2str(lamvec)]), xlabel('x'),ylabel('u(x)')
    subplot(1,2,2)
    plot(lamvec,mvec,'-*')
    title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
    drawnow, shg
end

% Set up initial tangent and tau we want to be orthogonal to.
told = chebfun(0,domain(uinit)); tauold = 1;

retract = 0; % retract == 1 if Newton told us to go back along the tangent.
if ( printing )
    if ( haveMeasure )
        fprintf('#Sol    #Newton     lambda     Steplength    Measure    \n')
        fprintf('--------------------------------------------------------\n')
    else
        fprintf('#Sol    #Newton     lambda     Steplength   \n')
        fprintf('--------------------------------------------\n')
    end
    %     fprintf('No. path iter    Newton iter   Steplength    Measure    Num. sol.\n')
    %     fprintf('----------------------------------------------------------------\n')
    
end
numSols = 1;
while counter <= maxiter
    % Find a tangent direction, but only if we were told by Newton not to
    % retract
    if ~retract
        [t, tau] = tangentBVP(N, uold, lamold, told, tauold, prefs);
        if counter == 1
            t = direction*t;
            tau = direction*tau;
        end
        % Move in the direction of the tangent
        uinit = uold+sl*t;
        laminit = lamold+sl*tau;
    end

    % Find a Newton correction
    [u, lam, iter, retract] = newtonBVP(N, uinit, laminit, t, tau, prefs);
    
    if retract % Newton told us we were trying to take too long tangent steps
        if printing 
            fprintf('retracted\n')
        end
        % Move in the direction of the current tangent, but only with
        % quarter of the steplength
        sl = sl/4;
        uinit = uold+sl*t;
        laminit = lamold+sl*tau;       
        if sl < smin
            disp('FAILED: sl < slmin')
            return
        end
        continue
    end
    
    if ( stopfun(u, lam) )
        return
    end
    % Store values for plotting
    
    if ( haveMeasure )
        measu = measure(u);
        mvec = [mvec; measu];
    end
    lamvec = [lamvec;lam];
    
    if printing
        if ( haveMeasure )
            fprintf('%3i \t   %2i \t   %6.2e \t %6.4f     %6.2e \n', ...
                counter,iter,lam, sl, measu)
        else
            fprintf('%3i \t   %2i \t   %6.2e \t %6.4f \n', ...
                counter,iter,lam, sl)
        end
    end
    
    % Count number of solutions found
%     if ( haveMeasure && mvec(end)*mvec(end-1) < 0 )
%         numSols = numSols + 1;
%     end


    if plotting
        subplot(1,2,1);
        plot(u),title(['Solution for \lambda =' num2str(lam)]), xlabel('x'),ylabel('u(x)')
        subplot(1,2,2)
        plot(lamvec,mvec,'-*'), 
        title('Bifurcation diagram'), xlabel('\lambda'), ylabel(mstring)
        drawnow, shg
    end
    
    counter = counter + 1;
    if iter >= 5
        sl = max(sl/2,smin); % Half steplength
    else
        sl = min(sl*2,smax); % Try to increase steplength
    end
    
    % If successful, update old values
    told = t; tauold = tau; uold = u; lamold = lam;
    
    % Update quasimatrix to be returned
    uquasi = [uquasi, u];
    
end


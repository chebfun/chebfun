function varargout = constructODEsol(solver, odefun, tspan, uinit, varargin)
%CONSTRUCTODESOL   Call one of Matlab's ODE solvers and return a CHEBFUN
%   Y = CHEBFUN.CONSTRUCTODESOL(SOLVER, ODEFUN, D, ...) calls one of the the
%   built-in MATLAB methods SOLVER to solve an initial-value problem on the
%   domain D. The result is then converted to a piecewise-defined CHEBFUN with
%   one column per solution component.
%
%   CHEBFUN.CONSTRUCTODESOL has the same calling sequence as Matlab's standard
%   odesolvers (ode113, ode15s, ode45) following the first argument.
%
%   One can also write [T, Y] = CONSTRUCTODESOL(...), in which case T is a
%   linear CHEBFUN on the domain D.
%
%   It is possible to pass a MATLAB ODESET struct to this method for specifying
%   options. The CHEBFUN overloads of the MATLAB ODE methods allow an extra
%   option, 'restartSolver', which if set to TRUE, will restart the ODE solver
%   at every breakpoint encountered. This is the default behaviour
%
% See also ODESET, ODE113, ODE15s, ODE45, ODE78, ODE89.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% By default, we want to restart the solver at breakpoints to catch behaviour
% that only happens at small intervals (cf.#1512). This can be overwritten by
% passing a 'restartSolver' option for the ODESET options structure passed to
% this method.
if ( (nargin > 4) && isfield(varargin{1}, 'restartSolver') )
    restartSolver = varargin{1}.restartSolver;
else
    restartSolver = true;
end


if ( ( length(tspan) == 2 ) || ~restartSolver )
    % We don't want to restart the solver, or we just have one interval.
    sol = solver(odefun, tspan, uinit, varargin{:});
    
    % Did we detect a blowup event?
    if ( isfield(sol,'ie') && ~isempty(sol.ie))
        oldEnd = tspan(2);
        tspan(2) = sol.xe; % Blowup time

        % Function for the solution where it's well defined:
        solFun = chebfun.odesol(sol, tspan, varargin{:});
        
        % Function of NaNs
        nanFun = chebfun(NaN(1, size(sol.y, 1)), [tspan(2), oldEnd]);

        % Joined fun:
        joinedFun = join(solFun, nanFun);
        
        % Did we request a time output as well?
        if ( nargout == 1 )
            varargout{1} = joinedFun;
        else
            timeFun = join(chebfun([tspan(1);tspan(2)],[tspan(1) tspan(2)]), ...
                chebfun([tspan(2); oldEnd], [tspan(2) oldEnd]));
            
            varargout{1} = timeFun;
            varargout{2} = joinedFun;
        end
    else
        [varargout{1:nargout}] = chebfun.odesol(sol, tspan, varargin{:});
    end
    
else
    % Here, we wish to restart the solver at each breakpoint encountered.
    
    % Number of pieces of the solution:
    numPieces = length(tspan) - 1;
    
    % Initialize a struct for storing the individual SOL pieces. Need extra
    % fields if event detection is on:
    %%% Birkisson suspects trouble here in older Matlab versions %%%
    if ( ~isempty(varargin) && ~isempty(varargin{1}.Events) )
        eventDetectionOn = true;
        sol(numPieces) = struct('solver', '', 'extdata', struct(), ...
            'x', [], 'y', [], 'stats', struct(), 'idata', struct(), ...
            'xe', [], 'ye', [], 'ie', []);
    else
        eventDetectionOn = false;
        sol(numPieces) = struct('solver', '', 'extdata', struct(), ...
            'x', [], 'y', [], 'stats', struct(), 'idata', struct());
    end
    
    % Keep track of whether the solver stopped due to event detection
    solverStopped = false;

    % Loop through the pieces
    for k = 1:numPieces
        % Compute the solution for the current interval
        sol(k) = solver(odefun, tspan(k:k+1), uinit, varargin{:});
        
        % Check if event detection is on, and if so, if we are to bail:
        if ( eventDetectionOn && ~isempty(sol(k).ie) )
            solverStopped = true;
            break
        else
            % Obtain a new initial condition for the next piece
            uinit = sol(k).y(:,end);
        end
    end
    
    % Convert all the pieces into a CHEBFUN:
    if ( ~solverStopped )
        % Easy case, event detection was off, or no stopping event detected:
        [varargout{1:nargout}] = chebfun.odesol(sol, tspan, varargin{:});
    else
        % Previously specified end time
        oldEnd = tspan(end);
        % Time periods where we have good solutions:
        tspan = tspan(1:k+1);
        
        % Blowup time
        tspan(end) = sol(k).xe; % Blowup time

        % Function for the solution where it's well defined:
        [timeFun, solFun] = chebfun.odesol(sol(1:k), tspan, varargin{:});
        
        % Function of NaNs
        nanFun = chebfun(NaN(1, size(sol(1).y, 1)), [tspan(end), oldEnd]);

        % Joined fun:
        joinedFun = join(solFun, nanFun);
        
        % Did we request a time output as well?
        if ( nargout == 1 )
            varargout{1} = joinedFun;
        else
            timeFun = join(timeFun, ...
                chebfun([tspan(2); oldEnd], [tspan(2) oldEnd]));
            
            varargout{1} = timeFun;
            varargout{2} = joinedFun;
        end
    end
end

end

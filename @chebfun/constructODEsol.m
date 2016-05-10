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
% See also ODESET, ODE113, ODE15s, ODE45.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


% By default, we want to restart the solver at breakpoints to catch behaviour
% that only happens at small intervals (cf.�#1512). This can be overwritten by
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
    [varargout{1:nargout}] = chebfun.odesol(sol, tspan, varargin{:});
    
else
    % Here, we wish to restart the solver at each breakpoint encountered.
    
    % Number of pieces of the solution:
    numPieces = length(tspan) - 1;
    
    % Initialize a struct for storing the individual SOL pieces:
    sol(numPieces) = struct('solver', '', 'extdata', struct(), ...
        'x', [], 'y', [], 'stats', struct(), 'idata', struct());

    % Loop through the pieces
    for k = 1:numPieces
        % Compute the solution for the current interval
        sol(k) = solver(odefun, tspan(k:k+1), uinit, varargin{:});
        % Obtain a new initial condition for the next piece
        uinit = sol(k).y(:,end);    
    end
    
    % Convert all the pieces into a CHEBFUN:  
    [varargout{1:nargout}] = chebfun.odesol(sol, tspan, varargin{:});
    
end

end

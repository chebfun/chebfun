function varargout = ode113(odefun, tspan, uinit, varargin)
%ODE113   Solve stiff differential equations and DAEs. Output a CHEBFUN.
%   Y = CHEBFUN.ODE113(ODEFUN, D, ...) applies the standard ODE113 method to
%   solve an initial-value problem on the domain D. The result is then converted
%   to a piecewise-defined CHEBFUN with one column per solution component.
%
%   CHEBFUN.ODE113 has the same calling sequence as Matlab's standard ODE113. 
%
%   One can also write [T, Y] = ODE113(...), in which case T is a linear CHEBFUN
%   on the domain D.
%
%   Note that CHEBFUN/ODE113() uses a default RELTOL of 1e-6.
%
%   It is possible to pass a MATLAB ODESET struct to this method for specifying
%   options. The CHEBFUN overloads of the MATLAB ODE methods allow an extra
%   option, 'resetSolver', which if set to TRUE, will restart the ODE solver at
%   every breakpoint encountered. This is the default behaviour
%
% Example:
%   y = chebfun.ode113(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:,1) - 1);                           % Find when y = 1
%
% See also ODESET, ODE15s, ODE45.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Option for disabling resetting 

% If varargin{1} (which correspond to options) does not have a field call
% resetSolver, or opts don't get passed in, this throws an error, which is
% resolved in the catch statement to go to the default behaviour of resetting.
try
    resetSolver = varargin{1}.resetSolver;
catch
    resetSolver = true;
end

% We don't want to reset the solver, or we just have one interval
if ( ( length(tspan) == 2 ) || ~resetSolver )
    sol = ode113(odefun, tspan, uinit, varargin{:});
    [t, y] = chebfun.odesol(sol, tspan, varargin{:});
else
    % Initialize a cell for storing the individual SOL pieces:
    solCell = cell(1, length(tspan) - 1);
    % Loop through the pieces
    for domCounter = 1:length(tspan)-1
        sol = ode113(odefun, tspan(domCounter:domCounter+1), uinit, varargin{:});
        solCell{domCounter} = sol;
        % Obtain a new initial condition for the next piece
        uinit = sol.y(:, end);
    end
    % Convert all the pieces into a CHEBFUN:
    [t, y] = chebfun.odesol(solCell, tspan, varargin{:});
end

% Output in a consistent way with the built in routine:
if ( nargout == 1 )
    % Only y will be returned in this case.
    varargout = {y};
else
    varargout = {t, y};
end

end

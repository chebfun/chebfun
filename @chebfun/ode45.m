function varargout = ode45(varargin)
%ODE45   Solve non-stiff differential equations. Output a CHEBFUN.
%   Y = CHEBFUN.ODE45(ODEFUN, D, ...) applies the standard ODE45 method to
%   solve an initial-value problem on the domain D. The result is then converted
%   to a piecewise-defined CHEBFUN with one column per solution component.
%
%   CHEBFUN.ODE45 has the same calling sequence as Matlab's standard ODE45. 
%
%   One can also write [T, Y] = ODE45(...), in which case T is a linear CHEBFUN
%   on the domain D.
%
%   Note that CHEBFUN/ODE45() uses a default RELTOL of 1e-6.
%
%   It is possible to pass a MATLAB ODESET struct to this method for specifying
%   options. The CHEBFUN overloads of the MATLAB ODE methods allow an extra
%   option, 'restartSolver', which if set to TRUE, will restart the ODE solver
%   at every breakpoint encountered.
%
% Example:
%   y = chebfun.ode45(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:, 1) - 1);                         % Find when y = 1
%
% See also ODESET, ODE113, ODE15S, ODE78, ODE89.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the CONSTRUCTODESOL method, with ode45 specified as the solver:
[varargout{1:nargout}] = chebfun.constructODEsol(@ode45, varargin{:});

end

function varargout = ode78(varargin)
%ODE78   Solve non-stiff differential equations. Output a CHEBFUN.
%   Y = CHEBFUN.ODE78(ODEFUN, D, ...) applies the standard ODE78 method to
%   solve an initial-value problem on the domain D. The result is then converted
%   to a piecewise-defined CHEBFUN with one column per solution component.
%
%   CHEBFUN.ODE78 has the same calling sequence as Matlab's standard ODE78.
%
%   One can also write [T, Y] = ODE78(...), in which case T is a linear CHEBFUN
%   on the domain D.
%
%   Note that CHEBFUN/ODE78() uses a default RELTOL of 1e-6.
%
%   It is possible to pass a MATLAB ODESET struct to this method for specifying
%   options. The CHEBFUN overloads of the MATLAB ODE methods allow an extra
%   option, 'restartSolver', which if set to TRUE, will restart the ODE solver
%   at every breakpoint encountered.
%
% Example:
%   y = chebfun.ode78(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:, 1) - 1);                         % Find when y = 1
%
% See also ODESET, ODE113, ODE15S, ODE45, ODE89.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the CONSTRUCTODESOL method, with ODE78 specified as the solver:
[varargout{1:nargout}] = chebfun.constructODEsol(@ode78, varargin{:});

end

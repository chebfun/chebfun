function varargout = ode15s(varargin)
%ODE15S   Solve stiff differential equations and DAEs. Output a CHEBFUN.
%   Y = CHEBFUN.ODE15S(ODEFUN, D, ...) applies the standard ODE15S method to
%   solve an initial-value problem on the domain D. The result is then converted
%   to a piecewise-defined CHEBFUN.
%
%   CHEBFUN.ODE15S has the same calling sequence as Matlab's standard ODE15S. 
%
%   One can also write [T, Y] = ODE15S(...), in which case T is a linear CHEBFUN
%   on the domain D.
%
%   Note that CHEBFUN/ODE15S() uses a default RELTOL of 1e-6.
%
%   It is possible to pass a MATLAB ODESET struct to this method for specifying
%   options. The CHEBFUN overloads of the MATLAB ODE methods allow an extra
%   option, 'restartSolver', which if set to TRUE, will restart the ODE solver
%   at every breakpoint encountered.
%
% Example:
%   y = chebfun.ode15s(@vdp1000, [0, 3000], [2; 0]); % Solve Van der Pol problem
%   roots(y(:,1) - 1);                               % Find when y = 1
%
% See also ODESET, ODE113, ODE45,

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the CONSTRUCTODESOL method, with ode15s specified as the solver:
[varargout{1:nargout}] = chebfun.constructODEsol(@ode15s, varargin{:});


end

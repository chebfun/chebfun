function varargout = ode113(varargin)
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
%   option, 'restartSolver', which if set to TRUE, will restart the ODE solver
%   at every breakpoint encountered. This is the default behaviour
%
% Example:
%   y = chebfun.ode113(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:,1) - 1);                           % Find when y = 1
%
% See also ODESET, ODE15s, ODE45.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the CONSTRUCTODESOL method, with ode113 specified as the solver:
[varargout{1:nargout}] = chebfun.constructODEsol(@ode113, varargin{:});

end

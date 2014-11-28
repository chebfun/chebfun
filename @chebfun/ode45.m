function varargout = ode45(varargin)
%ODE45   Solve stiff differential equations and DAEs. Output a CHEBFUN.
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
% Example:
%   y = chebfun.ode45(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:, 1) - 1);                         % Find when y = 1
%
% See also ODESET, ODE113, ODE15S,

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the built in ODE45():
sol = ode45(varargin{:});

% Convert solution to a CHEBFUN:
[t, y] = chebfun.odesol(sol, varargin{2}); 

% Output in a consistent way with the built in routine:
if ( nargout == 1 )
    % Only y will be returned in this case.
    varargout = {y};
else
    varargout = {t, y};
end

end

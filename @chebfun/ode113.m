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
% Example:
%   y = chebfun.ode113(@vdp1, [0, 20], [2 ; 0]); % Solve Van der Pol problem
%   roots(y(:,1) - 1);                           % Find when y = 1
%
% See also ODESET, ODE15s, ODE45.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Call the built in ODE113():
bigSol = struct('solver','ode113','x', [], 'y', []);
for domCounter = 1:length(tspan)-1
    sol = ode113(odefun, tspan(domCounter:domCounter+1), uinit, varargin{:});
    uinit = sol.y(:, end);
    % Convert solution to a CHEBFUN:
    [t, y] = chebfun.odesol(sol, tspan(domCounter:domCounter+1));
    ycell{domCounter} = y;
    tcell{domCounter} = t;
end

% Convert solution to a CHEBFUN:
% [t, y] = chebfun.odesol(sol, tspan); 
t = join(tcell{:});
y = join(ycell{:});
% Output in a consistent way with the built in routine:
if ( nargout == 1 )
    % Only y will be returned in this case.
    varargout = {y};
else
    varargout = {t, y};
end

end

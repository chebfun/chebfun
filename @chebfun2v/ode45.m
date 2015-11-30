function varargout = ode45(F,tspan,init,varargin)
% ODE45  Solve autonomous systems defined by a CHEBFUN2V.
%
%  [T, Y] = ODE45(F,TSPAN,Y0) with TSPAN = [T0 TFINAL] solves the autonomous
%  system of ODE y = f(y,y'), y'=g(y,y'), where f and g are the first and second
%  components of F, respectively, from time T0 to TFINAL with initial conditions
%  Y0. F is a CHEBFUN2V and Y is a complex-valued CHEBFUN representing the
%  solution, i.e., Y = y(t) + i*y'(t). To obtain solutions that interpolate at
%  T0,T1,...,TFINAL use TSPAN = [T0 T1 ... TFINAL].
%
%  [T,Y] = ODE45(F,TSPAN,Y0,OPTIONS) solves as above with default integration
%  properties replaced by values in OPTIONS, an argument created with the ODESET
%  function. See ODESET for details. However, the 'AbsTol' tolerance is always
%  set to machine precision.
%
%  [T,Y,TE,YE,IE] = ODE45(F,TSPAN,Y0,OPTIONS) with the 'Events' property in
%  OPTIONS set to a function handle EVENTS, solves as above while also finding
%  where functions of (T,Y), called event functions, are zero. For each function
%  you specify whether the integration is to terminate at a zero and whether the
%  direction of the zero crossing matters. These are the three column vectors
%  returned by EVENTS: [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y). For the I-th
%  event function: VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the
%  integration is to terminate at a zero of this event function and 0 otherwise.
%  DIRECTION(I)=0 if all zeros are to be computed (the default), +1 if only
%  zeros where the event function is increasing, and -1 if only zeros where the
%  event function is decreasing. Output TE is a column vector of times at which
%  events occur. Rows of YE are the corresponding solutions, and indices in
%  vector IE specify which event occurred.
%
%  SOL = ODE45(F,TSPAN,Y0,...) returns a structure storing information about
%  events. If events were detected, SOL.x is a row vector of points at which
%  events occurred. Columns of SOL.y are the corresponding solutions, and
%  indices in vector SOL.ie specify which event occurred.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO: When using events the default tolerance should be machine precision.]

if ( isempty( F ) )
    varargout = {[]};
    return
end

nF = F.nComponents;
dom = F.components{1}.domain; 

prefs = chebfunpref;
abstol = 100*prefs.techPrefs.eps;
% We don't expect 16 digits of relative tolerance, but at least require some.
reltol = 1e8*abstol;

% indicator function: 
ind = @(y) (y(1)>=dom(1)).*(y(1)<=dom(2)).*(y(2)>=dom(3)).*(y(2)<=dom(4));
g = @(t,y) feval(F, y(1), y(2)).*ind(y);
if ( nargin == 3 )
    opts = odeset('RelTol', reltol, 'AbsTol', abstol);
    sol = ode45(g, tspan, init, opts);
else
    opts = odeset('RelTol', reltol, 'AbsTol', abstol);
    opts = odeset(opts, varargin{:});
    sol = ode45(g, tspan, init, opts);
end


% Have a look at sol and see if events stopped the solution process. If so
% truncate the time interval.
tvec = sol.x;
if ( abs(tvec(end) - tspan(end)) > abstol )
    for jj = length(tspan):-1:1
        if tvec(end) <= tspan(jj)
            tspan(jj) = [];
        end
    end
    tspan = [tspan tvec(end)];
end

% Split solution at event locations.
if ( isfield(sol,'xe') )
    tvec = sol.xe;
    if ( ~isempty( tvec ) )
        tspan = [tspan tvec];
        tspan = unique( sort( tspan ) );
    end
end

t = chebfun(tspan([1 end]), tspan);
if ( any(any(isnan(sol.y))) )
    error('CHEBFUN:CHEBFUN2V:ode45:nan', ...
        'IVP returned NaN, try shorter time domain.')
else
    ys = sol.y;
end

if ( size(ys,2) == 2 )
    % Always ensure the result is complex-valued.
    y = chebfun(ys(:,1) + 1i*ys(:,2) + eps*1i, tspan);
else
    y = chebfun(ys(1,:).' + 1i*ys(2,:).' + eps*1i, tspan);
end

% We have no hope of storing this as a complex-valued CHEBFUN just pass back as
% a quasimatrix.
if ( nF == 3 )
    if size(ys,2) == 2
        y = chebfun(ys, tspan);
    else
        y = chebfun(ys.', tspan);
    end
end

if ( nargout > 1 )
    if ( nargout == 2 )
        varargout = {t, y};
    else
        varargout = {t, y, sol.x, sol.y, sol.ie};
    end
elseif ( nargout == 1 )
    cheb_sol = sol;
    cheb_sol.x = t;
    cheb_sol.y = y;
    varargout = {cheb_sol};
end

end

function [uout, tout] = spin2(varargin)
%SPIN2  Solve stiff PDEs in 2D periodic domains, Fourier spectral method and 
%exponential integrators.
%
%   UOUT = SPIN2(PDECHAR) solves the PDE specified by the string PDECHAR, and
%   plays a movie of the solution. Possible strings include 'GL' and 'GS' for 
%   the Ginzburg-Landau and Gray-Scott equations. Other PDEs are available, see 
%   Remark 1 and Examples 1-4. The output UOUT is a CHEBFUN2 corresponding to 
%   the solution at the final time (a CHEBMATRIX for systems of equations, each 
%   row representing one variable).
%
%   UOUT = SPIN2(S, N, DT) solves the PDE specified by the SPINOP2 S with N grid
%   points in each direction and time-step DT, and plays a movie of the solution. 
%   See HELP/SPINOP2 and Example 5.
%
%   [UOUT, TOUT] = SPIN2(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
%   Users of SPIN2 will quickly find they want to vary aspects of the plotting.
%   The fully general syntax for this involves using preferences specified by
%   a SPINPREF2 object PREF. See HELP/SPINPREF2 and Example 6. However for many 
%   purposes it is most convenient to use the syntax
%
%   UOUT = SPIN2(..., 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
%
%   For example:
%
%   UOUT = SPIN2(..., 'Clim', [a b]) changes colorbar limits to [a b] 
%   UOUT = SPIN2(..., 'colormap', 'jet') changes the colormap to 'jet'
%   UOUT = SPIN2(..., 'dataplot', 'abs') plots absolute value
%   UOUT = SPIN2(..., 'iterplot', 4) plots only every 4th time step 
%   UOUT = SPIN2(..., 'Nplot', 256) plays a movie at 256x256 resolution
%   UOUT = SPIN2(..., 'plot', 'off') for no movie
%   UOUT = SPIN2(..., 'view', [a b]) changes the view angle to [a b]
%
% Remark 1: List of PDEs (case-insensitive)
%
%    - 'GL' for the Ginzburg-Landau equation,
%    - 'GS' for the Gray-Scott equations (stripes),
%    - 'GSspots' for the Gray-Scott equations (spots),
%    - 'SCHNAK' for the Schnakenberg equations,
%    - 'SH' for the Swift-Hohenberg equation.
%
% Example 1: Ginzburg-Landau equation (spiral waves)
%
%       u = spin2('GL');
%
%    solves the Ginzburg-Landau equation
%
%       u_t = laplacian(u) + u - (1+1.5i)*u*|u|^2,
%
%    on [0 100]^2 from t=0 to t=100, with a RANDNFUN2 initial condition. 
%    The movie shows the real part of u.  For a movie of the absolute
%    value of u rather than the real part, execute
%
%       S = spinop2('GL');
%       u = spin2(S, 128, 1e-1, 'dataplot', 'abs');
%
% Example 2: Gray-Scott equations (pattern formation - stripes)
%
%       u = spin2('GS');
%
%    solves the Gray-Scott equations
%
%       u_t = 2e-5*laplacian(u) + F*(1-u) - u*v^2,
%       v_t = 1e-5*laplacian(v) - (F+K)*v + u*v^2,
%
%    on [0 1]^2 from t=0 to t=5000, with initial condition
%
%       u0(x,y) = 1 - exp(-100*((x-1/2.05)^2 + (y-1/2.05)^2)),
%       v0(x,y) = exp(-100*((x-1/2)^2 + 2*(y-1/2)^2)),
%
%    with F=0.030 and K=0.057.
% 
% Example 2 (bis): Gray-Scott equations (pattern formation - spots)
%
%       u = spin2('GSspots');
%
%    solves the Gray-Scott equations
%
%       u_t = 2e-5*laplacian(u) + F*(1-u) - u*v^2,
%       v_t = 1e-5*laplacian(v) - (F+K)*v + u*v^2,
%
%    on [0 1]^2 from t=0 to t=5000, with initial condition
%
%       u0(x,y) = 1 - exp(-100*((x-1/2.05)^2 + (y-1/2.05)^2)),
%       v0(x,y) = exp(-100*((x-1/2)^2 + 2*(y-1/2)^2)), 
%
%    with F=0.026 and K=0.059.
%
% Example 3: Schnakenberg equations (pattern formation - spots)
%
%       u = spin2('SCHNAK');
%
%    solves the Schnakenberg equations
%
%       u_t = laplacian(u) + 3*(a - u + u^2*v),
%       v_t = 10*laplacian(v) + 3*(b - u^2*v),
%
%    on [0 50]^2 from t=0 to t=500, with initial condition
%
%       u0(x,y) = (a+b) - exp(-2*((x-G/2.15)^2 + (y-G/2.15)^2)),
%       v0(x,y) = b/(a+b)^2 + exp(-2*((x-G/2)^2 + 2*(y-G/2)^2)),
%           with G=50, a=0.1 and b=0.9.
%
% Example 4: Swift-Hohenberg equation (Rayleigh-Benard convection rolls)
%
%       u = spin2('SH');
%
%    solves the Swift-Hohenberg equation
%
%       u_t = -2*laplacian(u) - biharmonic(u) - .9*u - u^3,
%
%    on [0 50]^2 from t=0 to t=800, with a RANDNFUN2 initial condition.
%
% Example 5: PDE specified by a SPINOP2
%
%       dom = [0 100 0 100]; tspan = [0 100];
%       S = spinop2(dom, tspan);
%       S.lin = @(u) lap(u);
%       S.nonlin = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
%       S.init = randnfun2(4, dom, 'trig');
%       S.init = S.init/norm(S.init, inf);
%       u = spin2(S, 128, 1e-1);
%
%   is equivalent to u = spin2('GL');
%
% Example 6: Using preferences
%
%       pref = spinpref2('plot', 'off', 'scheme', 'pecec433');
%       S = spinop2('SH');
%       u = spin2(S, 128, 5e-1, pref);
%   or simply,
%       u = spin2(S, 128, 5e-1, 'plot', 'off', 'scheme', 'pecec433');
%
%   solves the Swift-Hohenberg equation using N=128 grid points in each
%   direction, a time-step dt=5e-1, doesn't play any movie and uses the
%   time-stepping scheme PECEC433.
%
% See also SPINOP2, SPINPREF2, EXPINTEG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We are going to parse the inputs and call SOLVEPDE in the following ways,
%
%       SPINOPERATOR.SOLVEPDE(S, N, dt)
%  or
%       SPINOPERATOR.SOLVEPDE(S, N, dt, pref)
%
% where S is a SPINOP2 object, N is the number of grid points in each direction, 
% DT is the time-step and PREF is a SPINPREF2 oject.

% CASE 1. U = SPIN2('GL'):
if ( nargin == 1 ) 
    
    try spinop2(varargin{1});
    catch
        error('Unrecognized PDE. See HELP/SPINSPHERE for the list of PDEs.')
    end
    [S, N, dt, pref] = parseInputs(varargin{1});
    varargin{1} = S;
    varargin{2} = N;
    varargin{3} = dt;
    varargin{4} = pref;
    
% CASE 2. U = SPIN2('GL', 'PREF1', VALUE1) or U = SPINSPHERE(S, N, DT):
elseif ( nargin == 3 ) 
    
    % CASE 2.1. U = SPIN2('GL', 'PREF1', VALUE1):
    if ( isa(varargin{1}, 'char') == 1 && isa(varargin{2}, 'char') == 1 )
        [S, N, dt, pref] = parseInputs(varargin{1});
        pref.(varargin{2}) = varargin{3};
        varargin{1} = S;
        varargin{2} = N;
        varargin{3} = dt;
        varargin{4} = pref;
        
    % CASE 2.2. U = SPIN2(S, N, DT):
    else
        % Nothing to do here.
    end
    
% CASE 3. U = SPIN2(S, N, DT, PREF)
elseif ( nargin == 4 ) 
    % Nothing to do here.
    
% CASE 4. 
elseif ( nargin >= 5 )
    
    % CASE 4.1. U = SPIN2('GL', 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
    if ( isa(varargin{1}, 'char') == 1 && isa(varargin{2}, 'char') == 1 )
        [S, N, dt, pref] = parseInputs(varargin{1});
        j = 2;
        while j < nargin
            pref.(varargin{j}) = varargin{j+1};
            varargin{j} = [];
            varargin{j+1} = [];
            j = j + 2;
        end
        varargin{1} = S;
        varargin{2} = N;
        varargin{3} = dt;
        varargin{4} = pref;
        varargin = varargin(~cellfun(@isempty, varargin));
        
    % CASE 4.2. U = SPIN2(S, N, DT, 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
    else
        pref = spinpref2();
        j = 4;
        while j < nargin
            pref.(varargin{j}) = varargin{j+1};
            varargin{j} = [];
            varargin{j+1} = [];
            j = j + 2;
        end
        varargin{4} = pref;
        varargin = varargin(~cellfun(@isempty, varargin));
    end
    
end

% SPIN2 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end

function [S, N, dt, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse the inputs.

pref = spinpref2();
S = spinop2(pdechar);
if ( strcmpi(pdechar, 'GL') == 1 )
    dt = 1e-1;
    N = 128;
    pref.Clim = [-1 1];
    pref.iterplot = 2;
    pref.Nplot = 256;
elseif ( strcmpi(pdechar, 'GS') == 1 )
    dt = 5;
    N = 64;
    pref.Clim = [.25 .8 0 .35];
    pref.iterplot = 10;
    pref.Nplot = 128;
elseif ( strcmpi(pdechar, 'GSspots') == 1 )
    dt = 2;
    N = 64;
    pref.Clim = [.25 .85 0 .35];
    pref.iterplot = 25;
    pref.Nplot = 128;
elseif ( strcmpi(pdechar, 'SCHNAK') == 1 )
    dt = 5e-1;
    N = 64;
    pref.Clim = [.7 1.7 .65 1.05];
    pref.iterplot = 10;
    pref.Nplot = 128;
elseif ( strcmpi(pdechar, 'SH') == 1 )
    dt = 1;
    N = 128;
    pref.Clim = [-.4 .5];
    pref.iterplot = 4;
    pref.Nplot = 256;
end

end

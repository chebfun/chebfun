function [uout, tout] = spin3(varargin)
%SPIN3  Solve a time-dependent PDE in 3D with periodicity in space, using a
%Fourier spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN3(PDECHAR) solves the PDE specified by the string PDECHAR, and
%   plots a movie of the solution as it computes it; it is a demo mode.
%   The space and time intervals and the initial condition are chosen to produce
%   beautiful movies. Strings available include 'GL3' for Ginzburg-Landau
%   equation and 'GS3' for Gray-Scott equations. Many other PDEs are available,
%   see Remark 1 and Examples 1-4. The output UOUT is a CHEBFUN3 corresponding
%   to the solution at the final time (a CHEBMATRIX for systems of equations,
%   each row representing one variable).
%
%   UOUT = SPIN3(S, N, DT) solves the PDE specified by the SPINOP3 S with N grid
%   points in each direction and time-step DT. It plots a movie of the solution
%   as it computes it. See HELP/SPINOP3 and Example 5.
%
%   UOUT = SPIN3(S, N, DT, PREF) allows one to use the preferences specified by
%   the SPINPREF3 object PREF. See HELP/SPINPREF3 and Example 6.
%
%   UOUT = SPIN3(S, N, DT, 'PREF1', VALUEPREF1, 'PREF2', VALUEPREF2, ...) is an
%   alternative to the previous syntax. See Example 6.
%
%   [UOUT, TOUT] = SPIN3(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'GL3' for Ginzburg-Landau equation,
%    - 'GS3' for Gray-Scott equations,
%    - 'Schnak3' for Schnakenberg equations,
%    - 'SH3' for Swift-Hohenberg equation.
%
% Example 1: Ginzburg-Landau equation (spiral waves)
%
%       u = spin3('GL3');
%
%    solves the Ginzburg-Landau equation
%
%        u_t = laplacian(u) + u - (1+1.5i)*u*|u|^2,
%
%    on [0 50]^3 from t=0 to t=100, with a random initial condition.
%    The movie plots the real part of u.
%
% Example 2: Gray-Scott equations (pattern formation - fingerprints)
%
%       u = spin3('GS3');
%
%    solves the Gray-Scott equations
%
%       u_t = 3e-4*laplacian(u) + 3.5e-2*(1-u)*u - u*v^2,
%       v_t = 1.5e-4*laplacian(v) - 9.5e-2*v + u*v^2,
%
%    on [0 2]^3 from t=0 to t=6000, with initial condition
%
%       u0(x,y,z) = 1 - exp(-100*((x-G/2.05)^2 + (y-G/2.05)^2 + (z-G/2.15)^2)),
%       v0(x,y,z) = exp(-100*((x-G/2)^2 + 2*(y-G/2)^2 + 2*(z-G/2)^2)), with G=2.
%
% Example 3: Schnakenberg equations (pattern formation - spots)
%
%       u = spin3('Schnak3');
%
%    solves the Schnakenberg equations
%
%       u_t = laplacian(u) + 3*(.1 - u + u^2*v),
%       v_t = 10*laplacian(v) + 3*(.9 - u^2*v),
%
%    on [0 25]^3 from t=0 to t=500, with initial condition
%
%       u0(x,y,z) = (a+b) - exp(-2*((x-G/2.05)^2 + (y-G/2.05)^2 + (z-G/2.15)^2)),
%       v0(x,y,z) = b/(a+b)^2 + exp(-2*((x-G/2)^2 + 2*(y-G/2)^2 + 2*(z-G/2)^2)),
%           with G=25, a=0.1 and b=0.9.
%
% Example 4: Swift-Hohenberg equation (Rayleigh-Benard convection rolls)
%
%       u = spin3('SH3');
%
%    solves the Swift-Hohenberg equation
%
%       u_t = -2*laplacian(u) - biharmonic(u) - .9*u - u^3,
%
%    on [0 25]^3 from t=0 to t=800, with a random initial condition.
%
% Example 5: PDE specified by a SPINOP3
%
%       dom = [0 100 0 100 0 100]; tspan = [0 100];
%       S = spinop3(dom, tspan);
%       S.lin = @(u) lap(u);
%       S.nonlin = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
%       S.init = chebfun3(.1*randn(32, 32, 32), dom, 'trig')
%       u = spin3(S, 32, 1e-1);
%
%   is equivalent to u = spin3('GL3');
%
% Example 6: Using preferences
%
%       pref = spinpref3('plot', 'off', 'scheme', 'pecec433');
%       S = spinop3('sh3');
%       u = spin3(S, 32, 5e-1, pref);
%   or simply,
%       u = spin3(S, 32, 5e-1, 'plot', 'off', 'scheme', 'pecec433');
%
%   solves the Swift-Hohenberg equation using N=32 grid points in each
%   direction, a time-step dt=5e-1, doesn't produce any movie use the
%   time-stepping scheme PECEC433.
%
% See also SPINOP3, SPINPREF3, EXPINT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We are going to parse the inputs and call SOLVEPDE in the following ways,
%
%       SPINOPERATOR.SOLVEPDE(S, N, dt)
%  or
%       SPINOPERATOR.SOLVEPDE(S, N, dt, pref)
%
% where S is a SPINOP3 object, N is the number of grid points in each direction, 
% DT is the time-step and PREF is a SPINPREF3 oject.

if ( nargin == 1 ) % e.g., u = spin3('gl3')
    try spinop3(varargin{1});
    catch
        error('Unrecognized PDE. See HELP/SPIN3 for the list of PDEs.')
    end
    [S, N, dt, pref] = parseInputs(varargin{1});
    varargin{1} = S;
    varargin{2} = N;
    varargin{3} = dt;
    varargin{4} = pref;
elseif ( nargin == 3 ) % e.g., u = spin3(S, 32, 1e-1)
    % Nothing to do here.
elseif ( nargin == 4 ) % e.g., u = spin3(S, 32, 1e-1, pref)
    % Nothing to do here.
elseif ( nargin >= 5 ) % u.g., u = spin3(S, 32, 1e-1, 'plot', 'off')
    % In this case, put the options in a SPINPREF3 object.
    pref = spinpref3();
    j = 4;
    while j < nargin
        pref.(varargin{j}) = varargin{j+1};
        varargin{j} = [];
        varargin{j+1} = [];
        j = j + 2;
    end
    varargin{end + 1} = pref;
    varargin = varargin(~cellfun(@isempty, varargin));
end

% SPIN2 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end

function [S, N, dt, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse the inputs.

pref = spinpref3();
S = spinop3(pdechar);
if ( strcmpi(pdechar, 'GL3') == 1 )
    dt = 1e-1;
    N = 32;
    pref.Clim = [-1 1];
    pref.iterplot = 2;
    pref.Nplot = 64;
elseif ( strcmpi(pdechar, 'GS3') == 1 )
    dt = 5;
    N = 32;
    pref.Clim = [.3 .8 0 .35];
    pref.iterplot = 5;
    pref.Nplot = 64;
elseif ( strcmpi(pdechar, 'Schnak3') == 1 )
    dt = 5e-1;
    N = 32;
    pref.Clim = [.7 1.7 .65 1.05];
    pref.iterplot = 10;
    pref.Nplot = 64;
elseif ( strcmpi(pdechar, 'SH3') == 1 )
    dt = 1;
    N = 32;
    pref.Clim = [-.4 .5];
    pref.iterplot = 4;
    pref.Nplot = 64;
end

end
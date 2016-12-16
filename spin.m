function [uout, tout] = spin(varargin)
%SPIN  Solve a time-dependent PDE in 1D with periodicity in space, using a
%Fourier spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN(PDECHAR) solves the PDE specified by the string PDECHAR, and
%   plots a movie of the solution as it computes it; it is a demo mode.
%   The space and time intervals and the initial condition are chosen to produce
%   beautiful movies. Strings available include 'AC' for Allen-Cahn equation,
%   'KS' for Kuramoto-Sivashinsky equation and 'KdV' for Korteweg-de Vries
%   equation. Many other PDEs are available, see Remark 1 and Examples 1-9. The
%   output UOUT is a CHEBFUN corresponding to the solution at the final time
%   (a CHEBMATRIX for systems of equations, each row representing one variable).
%
%   UOUT = SPIN(S, N, DT) solves the PDE specified by the SPINOP S with N grid
%   points and a time-step DT. It plots a movie of the solution as it computes
%   it. See HELP/SPINOP and Example 10.
%
%   UOUT = SPIN(S, N, DT, PREF) allows one to use the preferences specified by
%   the SPINPREF object PREF. See HELP/SPINPREF and Example 11.
%
%   UOUT = SPIN(S, N, DT, 'PREF1', VALUEPREF1, 'PREF2', VALUEPREF2, ...) is an
%   alternative to the previous syntax. See Example 11.
%
%   [UOUT, TOUT] = SPIN(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'AC' for Allen-Cahn equation,
%    - 'Burg' for viscous Burgers equation,
%    - 'BZ' for Belousov-Zhabotinsky equation,
%    - 'CH' for Cahn-Hilliard equation,
%    - 'GS' for Gray-Scott equations,
%    - 'KdV' for Korteweg-de Vries equation,
%    - 'KS' for Kuramoto-Sivashinsky equation,
%    - 'Niko' for Nikolaevskiy equation,
%    - 'NLS' for the focusing nonlinear Schroedinger equation.
%
% Example 1: Allen-Cahn equation (metastable solutions)
%
%       u = spin('AC');
%
%    solves the Allen-Cahn equation
%
%       u_t = 5e-3*u_xx - u^3 + u,
%
%    on [0 2*pi] from t=0 to t=500, with initial condition
%
%       u0(x) = 1/3*tanh(2*sin(x)) - exp(-23.5*(x-pi/2)^2) + exp(-27*(x-4.2)^2)
%               + exp(-38*(x-5.4)^2).
%
% Example 2: Viscous Burgers equation (shock formation and dissipation)
%
%       u = spin('Burg');
%
%    solves the viscous Burgers equation
%
%           u_t = 1e-3*u_xx - u*u_x,
%
%   on [-1 1] from t=0 to t=15, with initial condition
%
%       u0(x) = (1-x^2)*exp(-30*(x+1/2)^2.
%
% Example 3: Belousov-Zhabotinsky (reaction-diffusion with three species)
%
%       u = spin('BZ');
%
%    solves the Belousov-Zhabotinsky equation
%
%       u_t = 1e-5*diff(u,2) + u + v - u*v - u^2,
%       v_t = 2e-5*diff(v,2) + w - v - u*v,
%       w_t = 1e-5*diff(w,2) + u - w
%
%    on [-1 1] from t=0 to t=30, with initial condition
%
%       u0(x) = exp(-100*(x+.5)^2),
%       v0(x) = exp(-100*(x)^2),
%       w0(x) = exp(-100*(x-.5)^2).
%
% Example 4: Cahn-Hilliard equation (metastable solutions)
%
%       u = spin('CH');
%
%    solves the Cahn-Hilliard equation
%
%       u_t = 1e-2*(-u_xx - 1e-3*u_xxxx + (u^3)_xx),
%
%    on [-1 1] from t=0 to t=100, with initial condition
%
%       u0(x) = 1/5*(sin(4*pi*x))^5 - 4/5*sin(pi*x).
%
% Example 5: Gray-Scott equations (pulse splitting)
%
%       u = spin('GS');
%
%    solves the Gray-Scott equations
%
%       u_t = diff(u,2) + 2e-2*(1-u) - u*v^2,
%       v_t = 1e-2*diff(u,2) - 8.62e-2*v + u*v^2,
%
%    on [-50 50] from t=0 to t=12000, with initial condition
%
%       u0(x) = 1 - 1/2*sin(pi*(x-L)/(2*L))^100,
%       v0(x) = 1/4*sin(pi*(x-L)/(2*L))^100,
%           with L=50.
%
% Example 6: Korteweg-de Vries equation (two-soliton solution)
%
%       u = spin('KdV');
%
%    solves the Korteweg-de Vries equation
%
%       u_t = -u*u_x - u_xxx,
%
%    on [-pi pi] from t=0 to t=0.03015, with initial condition
%
%       u0(x) = 3*A^2*sech(.5*A*(x+2))^2 + 3*B^2*sech(.5*B*(x+1))^2,
%           with A=25 and B=16.
%
% Example 7: Kuramoto-Sivashinsky (chaotic attractor)
%
%       u = spin('KS');
%
%    solves the Kuramoto-Sivashinsky equation
%
%       u_t = -u*u_x - u_xx - u_xxxx,
%
%    on [0 32*pi] from t=0 to t=300, with intial condition
%
%       u0(x) = cos(x/16)*(1 + sin(x/16)).
%
% Example 8: Nikolaevskiy equation (chaotic attractor)
%
%       u = spin('Niko');
%
%    solves the Nikolaevskiy equation
%
%       u_t = -u*u_x + 1e-1*u_xx + u_xxxx + u_xxxxxx,
%
%    on [0 32*pi] from t=0 to t=300, with intial condition
%
%       u0(x) = cos(x/16)*(1 + sin(x/16)).
%
% Example 9: Nonlinear Schroedinger equation (breather solution)
%
%       u = spin('NLS');
%
%    solves the focusing Nonlinear Schroedinger equation
%
%       u_t = i*u_xx + i*|u|^2*u,
%
%    on [-pi pi] from t=0 to t=18, with initial condition
%
%       u0(x) = 2*B^2/(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A,
%           with A=2 and B=1.%
%
%
% Example 10: PDE specified by a SPINOP
%
%       dom = [0 32*pi]; tspan = [0 300];
%       S = spinop(dom, tspan);
%       S.lin = @(u) -diff(u,2)-diff(u,4);
%       S.nonlin = @(u) -.5*diff(u.^2);
%       S.init = chebfun(@(x) cos(x/16).*(1 + sin(x/16)), dom, 'trig');
%       u = spin(S, 256, 5e-2);
%
%   is equivalent to u = spin('KS');
%
% Example 11: Using preferences
%
%       pref = spinpref('plot', 'waterfall', 'scheme', 'pecec433');
%       S = spinop('kdv');
%       u = spin(S, 256, 1e-5, pref);
%   or simply,
%       u = spin(S, 256, 1e-5, 'plot', 'waterfall', 'scheme', 'pecec433');
%
%   solves the KdV equation using N=256 grid points, a time-step dt=1e-5,
%   produces a WATERFALL plot as opposed to a movie and use the time-stepping
%   scheme PECEC433.
%
% See also SPINOP, SPINPREF, SPINSCHEME, SPIN2, SPIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We are going to parse the inputs and call SOLVEPDE in the following ways,
%
%       SPINOPERATOR.SOLVEPDE(S, N, dt)
%  or
%       SPINOPERATOR.SOLVEPDE(S, N, dt, pref)
%
% where S is a SPINOP object, N is the number of grid points, DT is the 
% time-step and PREF is a SPINPREF oject.

if ( nargin == 1 ) % e.g., u = spin('kdv')
    try spinop(varargin{1});
    catch
        error('Unrecognized PDE. See HELP/SPIN for the list of PDEs.')
    end
    [S, N, dt, pref] = parseInputs(varargin{1});
    varargin{1} = S;
    varargin{2} = N;
    varargin{3} = dt;
    varargin{4} = pref;
elseif ( nargin == 3 ) % e.g., u = spin(S, 256, 1e-5)
    % Nothing to do here.
elseif ( nargin == 4 ) % e.g., u = spin(S, 256, 1e-5, pref)
    % Nothing to do here.
elseif ( nargin >= 5 ) % u.g., u = spin(S, 256, 1e-5, 'plot', 'off')
    % In this case, put the options in a SPINPREF object.
    pref = spinpref();
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

% SPIN is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end

function [S, N, dt, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse the inputs.

pref = spinpref();
S = spinop(pdechar);
if ( strcmpi(pdechar, 'AC') == 1 )
    dt = 1e-1;
    N = 256;
elseif ( strcmpi(pdechar, 'Burg') == 1 )
    dt = 5e-3;
    N = 1024;
elseif ( strcmpi(pdechar, 'BZ') == 1 )
    dt = 1e-2;
    N = 256;
elseif ( strcmpi(pdechar, 'CH') == 1 )
    dt = 2e-2;
    N = 256;
elseif ( strcmpi(pdechar, 'GS') == 1 )
    dt = 5;
    N = 256;
elseif ( strcmpi(pdechar, 'KdV') == 1 )
    dt = 3e-6;
    N = 256;
elseif ( strcmpi(pdechar, 'KS') == 1 )
    dt = 5e-2;
    N = 256;
elseif ( strcmpi(pdechar, 'Niko') == 1 )
    dt = 1e-1;
    N = 256;
elseif ( strcmpi(pdechar, 'NLS') == 1 )
    dt = 1.5e-3;
    N = 256;
    pref.dataplot = 'abs';
end

end
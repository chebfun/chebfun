function [uout, tout] = spinsphere(varargin)
%SPINSPHERE  Solve a time-dependent PDE on the sphere with the double Fourier 
%sphere method and IMEX time-stepping algorithms.
%
%   UOUT = SPINSPHERE(PDECHAR) solves the PDE specified by the string PDECHAR,
%   and plots a movie of the solution as it computes it; it is a demo mode.
%   The space and time intervals and the initial condition are chosen to produce
%   beautiful movies. Strings available include 'AC' for Allen-Cahn equation, 
%   'GL' for Ginzburg-Landau equation, 'GM' for Gierer-Meinhardt equations and 
%   'NLS' for nonlinear Schrodinger equation. See Remark 1 and Examples 1-4. 
%   The output UOUT is a SPHEREFUN corresponding to the solution at the final 
%   time (a CHEBMATRIX for systems of equations, each row representing one 
%   variable).
%
%   UOUT = SPINSPHERE(S, N, DT) solves the PDE specified by the SPINOPSPHERE S 
%   with N grid points in each direction and time-step DT. It plots a movie of 
%   the solution as it computes it. See HELP/SPINOPSPHERE and Example 5.
%
%   UOUT = SPINSPHERE(S, N, DT, PREF) allows one to use the preferences 
%   specified by the SPINPREFSPHERE object PREF. See HELP/SPINPREFSPHERE and 
%   Example 6.
%
%   UOUT = SPINSPHERE(S, N, DT, 'PREF1', VALUEPREF1, 'PREF2', VALUEPREF2, ...) 
%   is an alternative to the previous syntax. See Example 6.
%
%   [UOUT, TOUT] = SPINSPHERE(...) also returns the times chunks TOUT at which 
%   UOUT was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'AC' for Allen-Cahn equation,
%    - 'GL' for Ginzburg-Landau equation,
%    - 'GM' for Gierer-Meinhardt equations,
%    - 'NLS' for focusing nonlinear Schroedinger equation.
%
% Example 1: Allen-Cahn equation
%
%        u = spinsphere('AC');
%
%    solves the Allen-Cahn equation
%
%        u_t = 1e-2*laplacian(u) + u - u^3
%
%    on the sphere from t=0 to t=60, with initial condition
%
%        u0(x, y, z) = cos(cosh(5*x.*z) - 10*y).
%
% Example 2: Ginzburg-Landau equation 
%
%        u = spinsphere('GL');
%
%    solves the Ginzburg-Landau equation
%
%        u_t = 5e-4*laplacian(u) + u - (1+1.5i)*u*|u|^2,
%
%    on the sphere from t=0 to t=100 with a RANDNFUNSPHERE initial condition.   
%    The movie plots the real part of u.
%
% Example 3: Gierer-Meinhardt equations
%
%        u = spinsphere('GM);
%
%    solves the Gierer-Meinhardt equations,
%
%       u_t = 1e-2*laplacian(u) + u.^2/v - u,
%       v_t = 1e-1*aplacian(v) + u^2 - v,
%
%    on the sphere from t=0 to t=40, with initial condition
%
%       u0(x,y,z) = v0(x,y,z) = 1+1/3*(cos(20*x)+cos(20*z)+cos(20*y)).
%
% Example 4: Focusing nonlinear Schroedinger equation
%
%        u = spinsphere('NLS');
%
%    solves the focusing nonlinear Schroedinger equation
%
%        u_t = 1i*laplacian(u) + 1i*u|u|^2,
%
%    on the sphere from t=0 to t=2, with initial condition
%
%     u0(lam, th) = (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*th)) - 1)*A 
%                  + 4*Y_6^6(lam, th), with A=2 and B=1.
%
%    The movie plots the real value of u.
%
% Example 4: PDE specified by a SPINOPSPHERE
%
%       tspan = [0 100];
%       S = spinopsphere(tspan);
%       S.lin = @(u) 1e-3*lap(u);
%       S.nonlin = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
%       S.init = spherefun(.1*randn(128));
%       u = spinsphere(S, 128, 2e-1);
%
%   is equivalent to u = spinsphere('GL');
%
% Example 5: Using preferences
%
%       pref = spinprefsphere('Clim', [-2 2]);
%       S = spinopsphere('ac');
%       u = spinsphere(S, 64, 2e-1, pref);
%   or simply,
%       u = spinsphere(S, 64, 2e-1, 'Clim', [-2 2]);
%
%   solves the Allen-Cahn equation using N=64 grid points in each direction,
%   a time-step dt=2e-1 and set the limits of the colobar to [-2 2].
%
% See also SPINOPSPHERE, SPINPREFSPHERE, IMEX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We are going to parse the inputs and call SOLVEPDE in the following ways,
%
%       SPINOPERATOR.SOLVEPDE(S, N, dt)
%  or
%       SPINOPERATOR.SOLVEPDE(S, N, dt, pref)
%
% where S is a SPINOPSPHERE object, N is the number of grid points in each 
% direction, DT is the time-step and PREF is a SPINPREFSPHERE oject.

if ( nargin == 1 ) % e.g., u = spinsphere('gl')
    try spinopsphere(varargin{1});
    catch
        error('Unrecognized PDE. See HELP/SPINSPHERE for the list of PDEs.')
    end
    [S, N, dt, pref] = parseInputs(varargin{1});
    varargin{1} = S;
    varargin{2} = N;
    varargin{3} = dt;
    varargin{4} = pref;
elseif ( nargin == 3 ) % e.g., u = spinsphere(S, 128, 1e-1)
    % Nothing to do here.
elseif ( nargin == 4 ) % e.g., u = spinsphere(S, 128, 1e-1, pref)
    % Nothing to do here.
elseif ( nargin >= 5 ) % u.g., u = spinsphere(S, 128, 1e-1, 'plot', 'off')
    % In this case, put the options in a SPINPREFSPHERE object.
    pref = spinprefsphere();
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

% SPINSPHERE is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end

function [S, N, dt, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse the inputs.

pref = spinprefsphere();
S = spinopsphere(pdechar);
if ( strcmpi(pdechar, 'AC') == 1 )
    dt = 1e-1;
    N = 256;
    pref.iterplot = 2;
elseif ( strcmpi(pdechar, 'GL') == 1 )
    dt = 1e-1;
    N = 256;
    pref.iterplot = 2;
elseif ( strcmpi(pdechar, 'GM') == 1 )
    dt = 1e-1;
    N = 128;
    pref.iterplot = 2;
elseif ( strcmpi(pdechar, 'NLS') == 1 )
    dt = 5e-3;
    N = 128;
    pref.iterplot = 2;
    pref.colormap = 'jet';
end

end
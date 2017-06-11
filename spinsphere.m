function [uout, tout] = spinsphere(varargin)
%SPINSPHERE  Solve stiff PDEs on the sphere, double Fourier sphere method and 
%implicit-explicit schemes.
%
%   UOUT = SPINSPHERE(PDECHAR) solves the PDE specified by the string PDECHAR,
%   and plays a movie of the solution. Possible strings include 'AC' and 'GL'
%   for the Allen-Cahn and Ginzburg-Landau equations. Other PDEs are available, 
%   see Remark 1 and Examples 1-4. The output UOUT is a SPHEREFUN corresponding
%   to the solution at the final time (a CHEBMATRIX for systems of equations, 
%   each row representing one variable).
%
%   UOUT = SPINSPHERE(S, N, DT) solves the PDE specified by the SPINOPSPHERE S 
%   with N grid points in each direction (longitude/latitude) and time-step DT, 
%   and plays a movie of the solution. See HELP/SPINOPSPHERE and Example 5.
%
%   UOUT = SPINSPHERE(S, N, DT, PREF) allows one to use the preferences 
%   specified by the SPINPREFSPHERE object PREF. See HELP/SPINPREFSPHERE and 
%   Example 6.
%
%   [UOUT, TOUT] = SPINSPHERE(...) also returns the times chunks TOUT at which 
%   UOUT was computed.
%
%   Users of SPINSPHERE will quickly find they want to vary aspects of the 
%   plotting. The fully general syntax for this involves using preferences 
%   specified by a SPINPREFSPHERE object PREF. See HELP/SPINPREFSPHERE and 
%   Example 6. However for many purposes it is most convenient to use the syntax
%
%   UOUT = SPINSPHERE(..., 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
%
%   For example:
%
%   UOUT = SPINSPHERE(..., 'Clim', [a b]) changes colorbar limits to [a b] 
%   UOUT = SPINSPHERE(..., 'colormap', 'jet') changes the colormap to 'jet'
%   UOUT = SPINSPHERE(..., 'dataplot', 'abs') plots absolute value
%   UOUT = SPINSPHERE(..., 'grid', 'on') for lagitude/longitude circles
%   UOUT = SPINSPHERE(..., 'iterplot', 4) plots only every 4th time step 
%   UOUT = SPINSPHERE(..., 'Nplot', 256) plays a movie at 256x256 resolution
%   UOUT = SPINSPHERE(..., 'plot', 'off') for no movie
%   UOUT = SPINSPHERE(..., 'view', [a b]) changes the view angle to [a b]
%
% Remark 1: List of PDEs (case-insensitive)
%
%    - 'AC' for the Allen-Cahn equation,
%    - 'GL' for the Ginzburg-Landau equation,
%    - 'GM' for the Gierer-Meinhardt equations,
%    - 'NLS' for the focusing nonlinear Schroedinger equation.
%
% Example 1: Allen-Cahn equation (metastable solutions)
%
%        u = spinsphere('AC');
%
%    solves the Allen-Cahn equation
%
%        u_t = 1e-2*laplacian(u) + u - u^3
%
%    on the sphere from t=0 to t=60, with initial condition
%
%        u0(x, y, z) = cos(cosh(5*x*z) - 10*y).
%
% Example 2: Ginzburg-Landau equation (spiral waves)
%
%        u = spinsphere('GL');
%
%    solves the Ginzburg-Landau equation
%
%        u_t = 1e-3*laplacian(u) + u - (1+1.5i)*u*|u|^2,
%
%    on the sphere from t=0 to t=100 with a RANDNFUNSPHERE initial condition.   
%    The movie shows the real part of u.
%
% Example 3: Gierer-Meinhardt equations (pattern formation - spots)
%
%        u = spinsphere('GM);
%
%    solves the Gierer-Meinhardt equations,
%
%       u_t = 1e-2*laplacian(u) + u^2/v - u,
%       v_t = 1e-1*laplacian(v) + u^2 - v,
%
%    on the sphere from t=0 to t=80, with initial condition
%
%       u0(x,y,z) = 1 + .1*(cos(20*x) + cos(20*y) + cos(20*z)),
%       v0(x,y,z) = 1 - .1*(cos(20*x) + cos(20*y) + cos(20*z)).
%
% Example 4: Nonlinear Schroedinger equation (spherical harmonic & breather)
%
%        u = spinsphere('NLS');
%
%    solves the focusing nonlinear Schroedinger equation
%
%        u_t = 1i*laplacian(u) + 1i*u|u|^2,
%
%    on the sphere from t=0 to t=3, with initial condition
%
%     u0(lam, th) = .1*(2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*th)) - 1)*A 
%                  + Y_8^6(lam, th), with A=1 and B=1.
%
%    The movie shows the absolute value of u.
%
% Example 5: PDE specified by a SPINOPSPHERE
%
%       tspan = [0 100];
%       S = spinopsphere(tspan);
%       S.lin = @(u) 1e-3*lap(u);
%       S.nonlin = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
%       S.init = randnfunsphere(.1);
%       S.init = S.init/norm(S.init, inf);
%       u = spinsphere(S, 128, 1e-1);
%
%   is equivalent to u = spinsphere('GL');
%
% Example 6: Using preferences
%
%       pref = spinprefsphere('Clim', [-1 1]);
%       S = spinopsphere('AC');
%       u = spinsphere(S, 128, 1e-1, pref);
%   or simply,
%       u = spinsphere(S, 128, 1e-1, 'Clim', [-1 1]);
%
%   solves the Allen-Cahn equation using N=128 grid points in each direction
%   and a time-step dt=1e-1, and sets the limits of the colorbar to [-1 1].
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
% direction, DT is the time-step and PREF is a SPINPREFSPHERE object.

% CASE 1. U = SPINSPHERE('GL'):
if ( nargin == 1 ) 
    
    try spinopsphere(varargin{1});
    catch
        error('Unrecognized PDE. See HELP/SPINSPHERE for the list of PDEs.')
    end
    [S, N, dt, pref] = parseInputs(varargin{1});
    varargin{1} = S;
    varargin{2} = N;
    varargin{3} = dt;
    varargin{4} = pref;
    
% CASE 2. U = SPINSPHERE('GL', 'PREF1', VALUE1) or U = SPINSPHERE(S, N, DT):
elseif ( nargin == 3 ) 
    
    % CASE 2.1. U = SPINSPHERE('GL', 'PREF1', VALUE1):
    if ( isa(varargin{1}, 'char') == 1 && isa(varargin{2}, 'char') == 1 )
        [S, N, dt, pref] = parseInputs(varargin{1});
        pref.(varargin{2}) = varargin{3};
        varargin{1} = S;
        varargin{2} = N;
        varargin{3} = dt;
        varargin{4} = pref;
        
    % CASE 2.2. U = SPINSPHERE(S, N, DT):
    else
        % Nothing to do here.
    end
    
% CASE 3. U = SPINSPHERE(S, N, DT, PREF)
elseif ( nargin == 4 ) 
    % Nothing to do here.
    
% CASE 4. 
elseif ( nargin >= 5 )
    
    % CASE 4.1. U = SPINSPHERE('GL', 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
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
        
    % CASE 4.2. U = SPINSPHERE(S, N, DT, 'PREF1', VALUE1, 'PREF2', VALUE2, ...)
    else
        pref = spinprefsphere();
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

% SPINSPHERE is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end

function [S, N, dt, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse the inputs.

pref = spinprefsphere();
S = spinopsphere(pdechar);
if ( strcmpi(pdechar, 'AC') == 1 )
    dt = 1e-1;
    N = 128;
    pref.Clim = [-1 1];
    pref.iterplot = 2;
    pref.Nplot = 256;
elseif ( strcmpi(pdechar, 'GL') == 1 )
    dt = 1e-1;
    N = 128;
    pref.Clim = [-1 1];
    pref.iterplot = 2;
    pref.Nplot = 256;
elseif ( strcmpi(pdechar, 'GM') == 1 )
    dt = 2e-1;
    N = 64;
    pref.Clim = [0 3 0.5 2];
    pref.iterplot = 4;
    pref.Nplot = 128;
elseif ( strcmpi(pdechar, 'NLS') == 1 )
    dt = 1e-2;
    N = 128;
    pref.colormap = 'jet';
    pref.dataplot = 'abs';
    pref.Clim = [0 1];
    pref.iterplot = 1;
    pref.Nplot = 256;
end

end

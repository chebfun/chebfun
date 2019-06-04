function varargout = pde15s(varargin)
%PDE15S   Solve PDEs using Chebfun.
%
%   UU = PDE15s(PDEFUN, T, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, T is a vector, U0 is a CHEBFUN or a CHEBMATRIX,
%   and BC is a CHEBOP boundary condition structure will solve the PDE dU/dt =
%   PDEFUN(U, t, x) with the initial condition U0 and boundary conditions BC
%   over the time interval T.
%
%   [TT, UU] = PDE15s(PDEFUN, T, U0, BC) returns also the time chunks TT.
%   If T is a 2-vector these will be given by linspace(T(1), T(2), 51).
%
%   [TT, UU, VV, ...] = PDE15s(...), in the case of a coupled system, will
%   return the quasimatrices UU, VV,... whose columns correspond to each
%   function at the times in the vector TT.
%
%   [TT, UU] = PDE15s(...), in the case of a coupled system, will return the
%   CHEBMATRIX UU, where each row corresponds to a variable at the times TT.
%
%   PDEFUN should take the form @(T, X, U1, U2, ..., UN), where U1, ..., UN are
%   the unknown dependent variables to be solved for, T is time, and X is space.
%
%   For backwards compatibility, the syntax @(U1, U2,..., UN, T, X, D, S, C)
%   for PDEFUN, where U1, ..., UN are the unknown dependent variables to be
%   solved for, T is time, X is space, D is the differential operator, S is the
%   definite integral operator (i.e., 'sum'), and C the indefinite integral
%   operator (i.e., 'cumsum'), is also supported.
%
%   By default, the solutions at intermediate time steps are not plotted when
%   PDE15s is called with output arguments. See documentation below on how
%   options can be specified and passed so that intermediate solutions get
%   plotted. If PDE15s is called without any output arguments, the solutions at
%   intermediate time steps are plotted by default.
%
%   For equations of one variable, UU is output as an array-valued CHEBFUN,
%   where UU(:, k) is the solution at T(k). For systems, the solution UU is
%   returned as a CHEBMATRIX with the different variables along the rows, and
%   time slices along the columns.
%
% Example 1: Nonuniform advection
%     x = chebfun('x', [-1 1]);
%     u = exp(3*sin(pi*x));
%     f = @(t, x, u) -(1 + 0.6*sin(pi*x)).*diff(u) + 5e-5*diff(u, 2);
%     opts = pdeset('Ylim', [0 20], 'PlotStyle', {'LineWidth', 2});
%     uu = pde23t(f, 0:.05:3, u, 'periodic', opts);
%     surf(uu, 0:.05:3)
%
% Example 2: Kuramoto-Sivashinsky
%     x = chebfun('x');
%     u = 1 + 0.5*exp(-40*x.^2);
%     bc.left = @(u) [u - 1 ; diff(u)];
%     bc.right = @(u) [u - 1 ; diff(u)];
%     f = @(u) u.*diff(u) - diff(u, 2) - 0.006*diff(u, 4);
%     opts = pdeset('Ylim', [-30 30], 'PlotStyle', {'LineWidth', 2});
%     uu = pde15s(f, 0:.01:.5, u, bc, opts);
%     surf(uu, 0:.01:.5)
%
% Example 3: Chemical reaction (system)
%      x = chebfun('x');
%      u = [ 1 - erf(10*(x+0.7)) ; 1 + erf(10*(x-0.7)) ; 0 ];
%      f = @(u, v, w)  [ .1*diff(u, 2) - 100*u.*v ; ...
%                        .2*diff(v, 2) - 100*u.*v ; ...
%                        .001*diff(w, 2) + 2*100*u.*v ];
%      opts = pdeset('Ylim', [0 2], 'PlotStyle', {'LineWidth', 2});
%      [t, u, v, w] = pde15s(f, 0:.1:3, u, 'neumann', opts);
%      mesh(uu{3})
%
% See chebfun/test/test_pde15s.m for more examples.
%
%   UU = PDE15s(PDEFUN, T, U0, BC, OPTS) will use nondefault options as defined
%   by the structure returned from OPTS = PDESET.
%
%   UU = PDE15s(PDEFUN, T, U0, BC, OPTS, N) will not adapt the grid size in
%   space. Alternatively OPTS.N can be set to the desired size.
%
%
%   There is some support for nonlinear and time-dependent boundary conditions,
%   such as
%       x = chebfun('x', [-1 1]);
%       u = exp(-3*x.^2);
%       f = @(t, x, u) .1*diff(u, 2);
%       bc.left = @(t, x, u) u - t;
%       bc.right = 0;
%       opts = pdeset('Ylim', [0 2], 'PlotStyle', {'LineWidth', 2});
%       uu = pde15s(f, 0:.1:2, u, bc, opts);
%       waterfall(uu);
%   with the input format being the same as PDEFUN described above.
%
% See also PDESET, ODE15S, PDE23T, PDESOLVE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Loop over the inputs and look for a PDESET struct:
optsPassed = false;                 % Set flag if we find optsd were passed.
for k = 1:nargin
    if ( isstruct(varargin{k}) && isfield(varargin{k}, 'ODESolver') )
        % Ammend the ODESolver field to be @ode15s:
        varargin{k}.ODESolver = @ode15s;
        optsPassed = true;          % We know now that options were passed.
        break                       % No need to continue the loop.
    end
end

% If no PDESET options were passed, create a PDESET struct:
if ( ~optsPassed )
    opts = pdeset('ODESolver', @ode15s);
    varargin{end+1} = opts;
end

% Call PDESOLVE() with option to use ODE15S():
[varargout{1:nargout}] = pdeSolve(varargin{:});

end

%% chebpde.m -- an executable m-file for solving a partial differential equation..
% Automatically created in CHEBGUI by user asgeir.
% Created on April 27, 2015 at 18:47.

%% Problem description.
% Solving
%   u_t = .03*u'' + u - u^3,
% for x in [-3,3] and t in [0,10], subject to periodic boundary conditions.

%% Problem set-up
% Create an interval of the space domain...
dom = [-3 3];
%...and specify a sampling of the time domain:
t = 0:0.2:10;

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) .03.*diff(u,2)+u-u.^3;

% Assign boundary conditions.
bc = 'periodic';

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial condition.
u0 = -1 + 2.*(exp(-35.*(x+2).^2) + exp(-11.*x.^2) + exp(-7.*(x-2).^2));

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'Plot', 'on','ylim',[-1 1]);

%% Call pde15s to solve the problem.
tic
[t, u] = pde15s(pdefun, t, u0, bc, opts);
toc
%% Plot the solution.
waterfall(u, t, 'simple', 'LineWidth', 2)
%% chebpdeSystem.m -- an executable m-file for solving a partial differential equation..
% Automatically created in CHEBGUI by user asgeir.
% Created on April 27, 2015 at 19:28.

%% Problem description.
% Solving
%   u_t = u" - v,
%   v" - u = 0,
% for x in [-1,1] and t in [0,2], subject to
%   u = 1, v = 1 at x = -1
% and
%   u = 1, v = 1 at x = 1

%% Problem set-up
% Create an interval of the space domain...
dom = [-1 1];
%...and specify a sampling of the time domain:
t = 0:.1:2;

% Make the right-hand side of the PDE.
% pdefun = @(u,v,t,x, diff) [diff(u,2)-v; -diff(v,2)+u];
pdefun = @(t,x, u,v) [diff(u,2)-v; -diff(v,2)+u];
pdeflag = [1  0]; % Zero when a variable is indep of time.

% Assign boundary conditions.
bc.left = @(t,u,v) [u-1; v-1];
bc.right = @(t,u,v) [u-1; v-1];

% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, dom);
% and of the initial conditions.
u0 = 1;
v0 = 1 - .5.*cos(.5.*pi.*x);
sol0 = [u0, v0];

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6, 'PDEflag', pdeflag, 'Ylim', [0.5,1], 'holdplot','on');

%% Call pde15s to solve the problem.
tic
[t, sol] = pde15s(pdefun, t, sol0, bc, opts);
toc

% Recover variable names.
u = sol(1,:);
v = sol(2,:);

%% Plot the solution.
figure
waterfall(sol, t, 'LineWidth', 2)
xlabel('x'), ylabel('t'), title('Solution')
figure
waterfall(u, t, 'LineWidth', 2)
xlabel('x'), ylabel('t'), title('u')
waterfall(v, t, 'LineWidth', 2)
xlabel('x'), ylabel('t'), title('v')
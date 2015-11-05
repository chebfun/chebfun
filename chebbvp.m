%% chebbvp.m -- an executable m-file for solving a boundary-value problem.
% Automatically created in CHEBGUI by user user.
% Created on August 09, 2015 at 02:33.

%% Problem description.
% Solving
%   u" - sin(v) = 0,
%   cos(u) + v" = 0,
% for x in [-1 1], subject to
%   u(-1)  = 1,
%   v'(-1) = 0,
%   u'(1)  = 0,
%   v(1)   = 0.

%% Problem set-up.

% Define the domain.
%Example 1:
%dom = [0 1];

%Example 2:
lambda = 1/3;
Bnd = 20;
%dom = [0 Bnd];
dom = [1 Bnd];

%K = 0.5; Re = 1e4; M =2;

sigma = 1; gamma = 3;

% Assign the differential equation to a chebop on that domain.
%N = chebop(@(x,u,v) [diff(u,2)-sin(v); cos(u)+diff(v,2)], dom);
%N = chebop(@(x,u,v) [(1+K)*diff(u,4)-Re*M*diff(u,2)+2*Re*u.*diff(u,3)-K*diff(v,2); (1+K/2)*diff(v,2)-Re*K*(2*v-diff(u,2))+Re*(2*u.*diff(v)-diff(u).*v)], dom);

%N = chebop(@(x,u,v) [diff(u,3) + v - diff(u).^2; diff(v,2)-3*sigma*diff(u).*v], dom);
N = chebop(@(x,u,v) [lambda^2*diff(u,2) + v - u.^2; lambda^2*diff(v,2)-3*sigma*u.*v], dom);

% Set up the rhs of the differential equation so that N(u, v) = rhs.
rhs = [0;0];

% Assign boundary conditions to the chebop.
%N.bc = @(x,u, v) [u(0); u(1); feval(diff(u),1)-1; feval(diff(u,2),0); v(1); v(0)];

%N.bc = @(x,u, v) [u(0); feval(diff(u),0); feval(diff(u),Bnd); v(0)-1; v(Bnd)];
N.bc = @(x,u, v) [u(1); u(Bnd); v(1)-1; v(Bnd)];

% Construct a linear chebfun on the domain, 
x = chebfun(@(x) x, dom);
% and assign an initial guess to the chebop.

% u_init = (x-1).^2 - 3;
% v_init = -(x+1).^2 + 4;
u_init = gamma*(1./x.^2 - 1./x.^3);
v_init = x.^-4;

N.init = [u_init, v_init];

%% Setup preferences for solving the problem.
% Create a CHEBOPPREF object for passing preferences.
% (See 'help cheboppref' for more possible options.)
options = cheboppref();

% Print information to the command window while solving:
options.display = 'iter';

% Option for tolerance.
options.errTol = 1e-16;
%options.errTol = 1e-3;

% Option for damping.
options.damping = false;

% Specify the discretization to use. Possible options are:
%  'collocation' (default)
%  'ultraspherical'
%  'periodic' (only for periodic problems).
options.discretization = 'collocation';
%options.discretization = 'ultraspherical';

% Option for determining how long each Newton step is shown.
options.plotting = 0.4;

%% Solve!
% Call solvebvp to solve the problem.
% (With the default options, this is equivalent to u = N\rhs.)
tic, [u, v] = solvebvp(N, rhs, options); toc

%% Plot the solution.
figure
%plot([u, v], 'LineWidth', 2)
%plot([diff(u)-1, v-1], 'LineWidth', 2)
%plot([u, v], 'LineWidth', 2)
%plot(x-1, u, x-1, v, 'LineWidth', 2)
plot(3*(x-1), u, 3*(x-1), v, 'LineWidth', 2)
% hold on
% plot(x-1, v, 'LineWidth', 2)
title('Final solution'), xlabel('x'), legend('u','v')
xlim([0 10]); shg

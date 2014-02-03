function [u1, u2, info1, info2] = test_scalarODE_breakpoints

% Setup
dom = [0 1 2 pi];
p = cheboppref;
p.plotting = 'on';
p.damped = 1;

N = chebop(@(x,u) diff(u,2) + sin(u-.2), dom);
N.lbc = @(u) u-2; N.rbc = @(u) u - 3;
N.init = chebfun(@(x) x/pi + 2, dom);
rhs = 0;

%% Try different discretizations
% Start with collocation -- no further action required
[u1, info1] = solvebvp(N, rhs, p);

%% Change to ultraS
p.discretization = @ultraS;
[u2, info2] = solvebvp(N, rhs, p);

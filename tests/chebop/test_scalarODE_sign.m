function [u1, u2, info1, info2] = test_scalarODE_sign

% Setup
dom = [-1 .5 1];
p = cheboppref;
p.plotting = 'on';
p.damped = 1;

N = chebop(@(x,u) diff(u,2) + sign(x).*sin(u), dom);
N.lbc = @(u) u - 2;
N.rbc = @(u) u - 2;
N.init = chebfun(@(x) 0*x + 2);
rhs = 0;

%% Try different discretizations
% Start with collocation -- no further action required
[u1, info1] = solvebvp(N, rhs, p);

%% Change to ultraS
p.discretization = @ultraS;
[u2, info2] = solvebvp(N, rhs, p);

function pass = test_firstOrderIntegralEqn(~)

%% Setup, testing a first order integral equation
d = [0 5];
L = chebop(@(u) diff(u) + 2*u + 5*cumsum(u), d);

%% Two different ways of specifying the LBC, should give the same results:
L.lbc = 0;
u1 = L\1;

L.lbc = [];
L.bc = @(x,u) u(0);
u2 = L\1;

%% Did we get the same solution?
pass = (norm(u1-u2) == 0);
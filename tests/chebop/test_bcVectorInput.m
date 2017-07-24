function pass = test_bcVectorInput(~)
% TEST_BCVECTORINPUT   Test specifying BCs with a vector

%% Nonlinear second order problem, first as an IVP
N = chebop(@(x,u) diff(u,2) + u.^2, [0.5 2]);
N.lbc = [0.1; -1];
u = N\0;
pass(1) = abs(u(.5) - 0.1) + abs(feval(diff(u), .5) + 1) < 5e-10;

%% Specify both lbc and rbc
N.lbc = .2;
N.rbc = -.2;
u = N\0;
pass(2) = abs(u(.5)-.2) + abs(u(2)+.2) < 1e-13;

%% Fourth order problem, first as an IVP
N = chebop(@(x,u) diff(u,4) + sin(u), [0.6 3]);
N.lbc = [0.1; 1; 3.2; -2];
u = N\0;
pass(3) = abs(u(.6) - 0.1) + abs(feval(diff(u), .6) - 1) + ...
    abs(feval(diff(u, 2), .6) - 3.2) + abs(feval(diff(u, 3), .6) + 2) < 2e-7;

%% Now as a BVP
N.lbc = [0.1 1];
N.rbc = [-.1 1];
u = N\0;
pass(4) = abs(u(.6) - 0.1) + abs(feval(diff(u), .6) - 1) + ...
    abs(u(3) + 0.1) + abs(feval(diff(u),3) - 1) < 1e-11;

%% Both in one go:
N.lbc = [];
N.rbc = [];
N.bc = [0.2 0.8];
u = N\0;
pass(5) = abs(u(.6) - 0.2) + abs(feval(diff(u), .6) - 0.8) + ...
    abs(u(3) - 0.2) + abs(feval(diff(u),3) - 0.8) < 1e-11;

%% First order systems are supported
N = chebop(@(x,u,v) [diff(u) + v; diff(v) - u]);
N.lbc = [1; 2];
[u,v] = N\0;
pass(6) = abs(u(-1) - 1) + abs(v(-1) - 2) < 1e-11;

%% First order system (more variables)
N = chebop(@(x,u,v,w,y) [diff(u) + v; diff(v) - w; diff(w)-y; diff(y)-u]);
N.rbc = [1; 2; 4; 3];
[u,v,w,y] = N\0;
pass(7) = abs(u(1) - 1) + abs(v(1) - 2) + abs(w(1) - 4) + ...
    abs(y(1) - 3) < 1e-11;

%% First order system (more variables), test rbc
N = chebop(@(x,u,v,w,y) [diff(u) + v; diff(v) - w; diff(w)-y; diff(y)-u]);
N.rbc = [1; 2; 4; 3];
[u,v,w,y] = N\0;
pass(8) = abs(u(1) - 1) + abs(v(1) - 2) + abs(w(1) - 4) + ...
    abs(y(1) - 3) < 1e-11;

%% Higher order systems are not supported
N = chebop(@(x,u,v) [diff(u,2) + v; diff(v) - u]);
N.lbc = [1; 3; 2];
try
    [u,v] = N\0;
catch ME
    pass(9) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:parseBC:numberOfConditions');
end
end
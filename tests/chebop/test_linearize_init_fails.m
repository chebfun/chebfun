function pass = test_linearize_init_fails(~)
%TEST_LINEARIZE_INIT_FAILS   Test that we get sensible error messages when we
%are dealing with an initial guess that makes the chebop fail to evaluate.

%% First order problem which is OK, as it's solved as an IVP
N = chebop(0, 1);
N.op = @(u) diff(u) - sqrt(u);
N.lbc = 1;
u = N\1;
pass(1) = norm(N(u)-1) < 1e-10;

%% First order problem which is OK, not OK when try to solve as a BVP
N = chebop(0, 1);
N.op = @(u) diff(u) - sqrt(u);
N.bc = @(u) u(0) - 1;
try
    u = N\1;
catch ME
    pass(2) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:linearize:invalidInitialGuess');
end

%% Nor when solvebvp is called directly
N = chebop(0, 1);
N.op = @(u) diff(u) - sqrt(u);
N.lbc = 1;
try
    u = solvebvp(N, 1);
catch ME
    pass(3) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:linearize:invalidInitialGuess');
end

%% Second order problem
N = chebop(0, 1);
N.op = @(u) diff(u,2) - sqrt(u);
N.bc = 1;
try
    u = N\1;
catch ME
    pass(4) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:linearize:invalidInitialGuess');
end

%% Second order problem, different function causing issues
N = chebop(0, 1);
N.op = @(u) diff(u,2) - 1./u;
N.bc = 1;
try
    u = N\1;
catch ME
    pass(5) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:linearize:invalidInitialGuess');
end

%% Coupled system, OK when solved as an IVP
N = chebop(0, 1);
N.op = @(x, u, v) [diff(u,2) - sqrt(v); diff(v,2) + 1./u];
N.lbc = @(u,v) [u-1; diff(u)-.1; v-2; diff(v) + .2];
uv = N\[1; 2];
pass(6) = norm(N(uv) - [1;2]) < 1e-8;

%% Coupled system, OK when solved as an IVP
N = chebop(0, 1);
N.op = @(x, u, v) [diff(u,2) - sqrt(v); diff(v,2) + 1./u];
N.bc = @(x,u,v) [u(0) - 1; u(1) + 1; v(0) - 2; v(0) - 3];
try
    uv = N\1;
catch ME
    pass(7) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBOP:linearize:invalidInitialGuess');
end
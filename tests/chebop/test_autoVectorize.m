function pass = test_autoVectorize
%TEST_AUTOVECTORIZE    Test the automatic vectorization works for CHEBOPs
%
% Note: This test only checks whether things run, not the correctness of the
% results.
%% Setup
% Domain and some of CHEBFUNs
d = [0, 1];
x = chebfun(@(x) x, d);
f = exp(sin(4*x));
g = cos(x.^2).^2;
h = 1./(5+x.^2);
u = sin(5*x.^2);
v = exp(x/2);
% Some parameters
a = 5;
b = 4*pi;
c = -5*exp(1);

%% Scalar problem -- IVP, pass OP to constructor.
fun = @(u) diff(u, 2) - a*(1-u^2)*diff(u) + u;
N = chebop(fun, d);
N.lbc = @(u) [u - 2; diff(u)];
try
    u = N\0;
    pass(1) = 1;
catch ME
    pass(1) = 0;
end

%% Scalar problem -- IVP, pass OP after constructing.
fun = @(x,u) diff(u, 2) - f*(1-u^2)*diff(u) + u;
N = chebop(d);
N.op = fun;
N.lbc = @(u) [u - 2; diff(u)];
try
    u = N\0;
    pass(2) = 1;
catch ME
    pass(2) = 0;
end

%% Scalar problem -- BVP, pass OP to constructor.
fun = @(u) diff(u, 2) - a*(1-u^2)*diff(u) + u;
N = chebop(fun, d);
N.bc = @(x,u) [u(0) - 2; u(1) - 3];
try
    u = N\0;
    pass(3) = 1;
catch ME
    pass(3) = 0;
end

%% Scalar problem -- BVP, pass OP after constructing.
fun = @(x,u) diff(u, 2) - .2*f*(1-u^2)*diff(u) + u;
N = chebop(d);
N.op = fun;
N.bc = @(x,u) [u(0) - 2; u(1) - 3];
try
    u = N\0;
    pass(4) = 1;
catch ME
    pass(4) = 0;
end

%% Coupled problem -- BVP, pass OP to constructor, also need to vectorize BCs
fun = @(x,u,v) [diff(u, 2) - a*u^2*diff(v); diff(v,2) + diff(u)/(v+2)^2];
N = chebop(fun, d);
N.bc = @(x,u,v) [u(0) - 2; u(1) - 3; v(0)-1; v(1)*feval(diff(u), 1) - 17.909];
u0 = chebfun([2.0000 1.5059  1.2210 2.0354 3.0000]', d);
v0 = chebfun([1; 2], d);
N.init = [u0; v0];
try
    uv = N\0;
    pass(5) = 1;
catch ME
    pass(5) = 0;
end

%% Scalar problem -- BVP. Vectorization off, should give an error
fun = @(u) diff(u, 2) - u*diff(u);
N = chebop(d);
N.vectorize = false;
N.op = fun;
N.bc = @(x,u) [u(0) - 2; u(1) - 3];
try
    u = N\0;
    % We should have been able to solve this!
    pass(6) = 0;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:ADCHEBFUN:mtimes:dims');
end

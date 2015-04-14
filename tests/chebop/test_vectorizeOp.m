function pass = test_vectorizeOp()
%TEST_VECTORIZEOP    Test that the CHEBOP vectorization method works as expected

%% Setup
% Domain and some of CHEBFUNs
d = [1, 3];
x = chebfun(@(x) x, d);
f = exp(sin(4*x));
g = cos(x.^2).^2;
h = 1./(5+x.^2);
u = sin(5*x.^2);
v = exp(x/2);
% Some parameters
b = 4*pi;
c = -5*exp(1);

%% Simple case, already vectorized
fun = @(u) diff(u, 2) + b.*u.^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect (FUN above is already vectorized, so they're identical)
exactFun = fun;
% Check that we got the correct string and the correct function
pass(1,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(1,2) = norm(exactFun(u) - vecFun(u)) == 0;

%% Simple case, already vectorized, with both x and u
fun = @(x,u) diff(u, 2) + b.*u.^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect (FUN above is already vectorized, so they're identical)
exactFun = fun;
% Check that we got the correct string and the correct function
pass(2,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(2,2) = norm(exactFun(x,u) - vecFun(x,u)) == 0;

%% Vectorization needed, scalar case
fun = @(u) diff(u, 2) + b*u^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(u) diff(u, 2) + b.*u.^3;
% Check that we got the correct string and the correct function
pass(3,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(3,2) = norm(exactFun(u) - vecFun(u)) == 0;

%% Vectorization needed, scalar case with both x and u
fun = @(x,u) diff(u, 2) + b*u^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(x,u) diff(u, 2) + b.*u.^3;
% Check that we got the correct string and the correct function
pass(4,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(4,2) = norm(exactFun(x,u) - vecFun(x,u)) == 0;

%% Vectorization needed, scalar case with a CHEBFUN
fun = @(u) diff(u, 2) + f*u^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(u) diff(u, 2) + f.*u.^3;
% Check that we got the correct string and the correct function
pass(5,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(5,2) = norm(exactFun(u) - vecFun(u)) == 0;

%% Vectorization needed, scalar case with a CHEBFUN, both x and u
fun = @(x,u) diff(u, 2) + f*u^3;
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(x,u) diff(u, 2) + f.*u.^3;
% Check that we got the correct string and the correct function
pass(6,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(6,2) = norm(exactFun(x,u) - vecFun(x,u)) == 0;

%% System case, already vectorized
fun = @(x,u,v) [diff(u, 2) + b.*u.^3+v.*u./h; diff(v,2) + v.^2./c + g.*u./v];
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect (FUN above is already vectorized, so they're identical)
exactFun = fun;
% Check that we got the correct string and the correct function
pass(7,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(7,2) = norm(exactFun(x,u,v) - vecFun(x,u,v)) == 0;

%% System case, partially vectorized
fun = @(x,u,v) [diff(u, 2) + b*u.^3+v.*u/h; diff(v,2) + v.^2./c + g*u./v];
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(x,u,v) [diff(u, 2) + b.*u.^3+v.*u./h; ...
    diff(v,2) + v.^2./c + g.*u./v];
% Check that we got the correct string
pass(8,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(8,2) = norm(exactFun(x,u,v) - vecFun(x,u,v)) == 0;

%% System case, no vectorization
fun = @(x,u,v) [diff(u, 2) + b*u^3+v*u/h; diff(v,2) + v^2/c + g*u/v];
% Automatically vectorized
vecFun = chebop.vectorizeOp(fun);
% The answer we expect
exactFun = @(x,u,v) [diff(u, 2) + b.*u.^3+v.*u./h; ...
    diff(v,2) + v.^2./c + g.*u./v];
% Check that we got the correct string and the correct function
pass(9,1) = strcmp(func2str(exactFun), func2str(vecFun));
pass(9,2) = norm(exactFun(x,u,v) - vecFun(x,u,v)) == 0;

%% For a function handle, will get the same results
s = @sin;
sVec = chebop.vectorizeOp(s);
% Check that we got the correct string and the correct function
pass(10,1) = strcmp(func2str(s), func2str(sVec));
pass(10,2) = norm(s(u)-sVec(u)) == 0;
end
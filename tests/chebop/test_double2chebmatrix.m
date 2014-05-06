% TEST_DOUBLE2CHEBMATRIX
function pass = test_double2chebmatrix
% Tests passing various arguments to CHEBOP/double2chebmatrix()
% AB, 05/05/2014

% Note that below, we check whether the error in the output from
% double2chebmatrix() is precisely 0. That is fine, since that method should
% return CHEBMATRIX objects identical to what we compare them with.

%% Scalar case
% Setup
dom = [0, 2];
x = chebfun(@(x) x, dom);

% Dummy residual function
res = chebmatrix({x});

out = chebop.double2chebmatrix(4, res);

pass(1) = ( norm(out{1} - 4) == 0 );

%% System case, two CHEBFUNS
res = [x; x];
out = chebop.double2chebmatrix([3; 5], res);
pass(2) = ( norm(out - [0*x + 3; 0*x + 5]) == 0 );

%% System case, CHEBFUN and a double
res = [x; 2];
out = chebop.double2chebmatrix([3; 5], res);
pass(3) = ( norm(out - [0*x + 3; 5]) == 0 );

%% Same as above, but now with a breakpoint!
dom = [0, .6, 2];
x = chebfun(@(x) x, dom);

% Dummy residual function
res = chebmatrix({x});

out = chebop.double2chebmatrix(4, res);

pass(4) = ( norm(out{1} - 4) == 0 );

%% System case, two CHEBFUNS
res = [x; x];
out = chebop.double2chebmatrix([3; 5], res);
pass(5) = ( norm(out - [0*x + 3; 0*x + 5]) == 0 );

%% System case, CHEBFUN and a double
res = [x; 2];
out = chebop.double2chebmatrix([3; 5], res);
pass(6) = ( norm(out - [0*x + 3; 5]) == 0 );
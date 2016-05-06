% Test file for chebtech/clenshaw.m

function pass = test_clenshaw(varargin)

% Set a tolerance (pref.chebfuneps doesn't matter)
tol = 10*eps;

%%
% Test that a single coefficient is evaluated correctly:
% For a scalar evaluation:
c = sqrt(2);
v = chebtech.clenshaw(0, c);
pass(1) = c == v;

% For a vector evaluation:
x = [-.5 ; 1];
v = chebtech.clenshaw(x, c);
pass(2) = ( all(size(v) == [2, 1]) && all(c == v) );

% For a row vector evaluation with column coefficients:
c = [c, c, c];
v = chebtech.clenshaw(x, c);
pass(3) = ( all(size(v) == [2, 3]) && norm(repmat(c, 2, 1) - v) == 0);

%%
% Test that a vector coefficient is evaluated correctly:
% Some simple data :
c = (5:-1:1).';
x = [-.5 ; -.1 ; 1];

% Scalar coefficient
v = chebtech.clenshaw(x, c);
% Exact values:
vTrue = [3 ; 3.1728 ; 15];
pass(4) = norm(v - vTrue, inf) < tol;

% In vectorised form:
vTrue2 = [0 ; 3.6480 ; 15];
v = chebtech.clenshaw(x, [c, c(end:-1:1)]);
pass(5) = norm(v - [vTrue, vTrue2], inf) < tol;

end

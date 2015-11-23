function pass = test_deriv(~)
% TEST_DERIV   Test the chebfun/deriv method

% This test creates some chebfun, evaluates the derivative and compares with the
% results of feval(diff(...)).

% Our CHEBFUN
u = chebfun(@(x) sin(exp(x)), [0 2]);

%% Evaluation at one point, default value of derivative
pass(1) = norm(deriv(u, 1) - feval(diff(u), 1)) == 0;

%% Higher order derivative
pass(2) = norm(deriv(u, .5, 3) - feval(diff(u, 3), .5)) == 0;

%% Default derivative at a vector of points
xx = linspace(0.1, 0.5, 11);
pass(3) = norm(deriv(u, xx) - feval(diff(u), xx)) == 0;

%% Higher derivative, vector of points
pass(4) = norm(deriv(u, xx, 4) - feval(diff(u, 4), xx)) == 0;

%% Array valued CHEBFUN at a vector of points
pass(5) = norm(deriv([u cos(u)], xx) - feval(diff([u cos(u)]), xx)) == 0;
pass(6) = norm(deriv([u cos(u)], xx, 2) - feval(diff([u cos(u)], 2), xx)) == 0;

end
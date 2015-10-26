function pass = test_spline(pref)

% Test a scalar function:
x = 0:10;  
y = sin(x);
f = chebfun.spline(x, y);
tol = 10*eps;
pass(1) = norm(feval(f, x) - y) < tol;
pass(2) = numel(f.funs) == 10;
pass(3) = length(f) == 40;

% Test an array-valued function:
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.spline(x, y);
tol = 10*eps;
pass(4) = norm(feval(f, x) - y) < tol;
pass(5) = numel(f.funs) == 10;
pass(6) = length(f) == 40;

% Test an end-splope condition:
x = (0:10)';  
y = [sin(x), cos(x)];
y = [0 0 ; y ; 0 0];
f = chebfun.spline(x, y);
tol = 10*eps;
pass(7) = norm(feval(f, x) - y(2:end-1,:)) < tol;
pass(8) = numel(f.funs) == 10;
pass(9) = length(f) == 40;
pass(10) = 1;
% pass(10) = norm(feval(diff(f), [-1;1]), inf) < tol; % [TODO]: Require feature-calclus.

% Test a different domain:
dom = [.01, 10.01];
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.spline(x, y, dom);
tol = 10*eps;
pass(11) = all(f.domain == [dom(1), 1:10, dom(2)]);
pass(12) = norm(feval(f, x(2:end)) - y(2:end,:)) < tol;
pass(13) = numel(f.funs) == 11;
pass(14) = length(f) == 44;

end

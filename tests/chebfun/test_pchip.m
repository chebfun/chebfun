function pass = test_pchip(pref)

% Test a scalar function:
x = 0:10;  
y = sin(x);
f = chebfun.pchip(x, y);
tol = 10*eps;
pass(1) = norm(feval(f, x) - y) < tol;
pass(2) = numel(f.funs) == 10;
pass(3) = length(f) == 40;

% Test an array-valued function:
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.pchip(x, y);
tol = 10*eps;
pass(4) = norm(feval(f, x) - y) < tol;
pass(5) = numel(f.funs) == 10;
pass(6) = length(f) == 40;

% Test a different domain:
dom = [.01, 10.01];
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.pchip(x, y, dom);
tol = 10*eps;
pass(7) = all(f.domain == [dom(1), 1:10, dom(2)]);
pass(8) = norm(feval(f, x(2:end)) - y(2:end,:)) < tol;
pass(9) = numel(f.funs) == 11;
pass(10) = length(f) == 44;

% Test the example from the m-file:
x = -3:3;
y = [-1 -1 -1 0 1 1 1];
f = chebfun.pchip(x, y);
tol = 10*eps;
pass(11) = norm(feval(f, x) - y) < tol;
pass(12) = numel(f.funs) == 6;

end

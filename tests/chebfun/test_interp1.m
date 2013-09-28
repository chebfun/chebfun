function pass = test_interp1(pref)

%% Linear interpolation:
% Test a scalar function:
x = (0:10).';  
y = sin(x);
f = chebfun.interp1(x, y, 'linear');
tol = epslevel(f);
pass(1) = norm(feval(f, x) - y) < tol;
pass(2) = numel(f.funs) == 10;
pass(3) = length(f) == 20;

% Test an array-valued function:
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.interp1(x, y, 'linear');
tol = epslevel(f);
pass(4) = norm(feval(f, x) - y) < tol;
pass(5) = numel(f.funs) == 10;
pass(6) = length(f) == 20;

% Test a different domain:
dom = [.01, 10.01];
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.interp1(x, y, 'linear', dom);
tol = epslevel(f);
pass(7) = all(f.domain == [dom(1), 1:10, dom(2)]);
pass(8) = norm(feval(f, x(2:end)) - y(2:end,:)) < tol;
pass(9) = numel(f.funs) == 11;
pass(10) = length(f) == 22;

%% Polynomial interpolation:
% Test an scalar-valued function:
x = (0:10)';  
y = sin(x);
f = chebfun.interp1(x, y);
tol = 100*epslevel(f);
pass(11) = norm(feval(f, x) - y) < tol;
pass(12) = numel(f.funs) == 1;
pass(13) = length(f) == length(x);

% Test an array-valued function:
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.interp1(x, y);
tol = 100*epslevel(f);
pass(14) = norm(feval(f, x) - y) < tol;
pass(15) = numel(f.funs) == 1;
pass(16) = length(f) == length(x);

% Test a different domain:
dom = [.01, 10.01];
f = chebfun.interp1(x, y, dom);
tol = 100*epslevel(f);
pass(17) = norm(feval(f, x) - y) < tol;
pass(18) = numel(f.funs) == 1;
pass(19) = length(f) == length(x);
pass(20) = all(f.domain == dom);

end

function pass = test_interp1(pref)

%% Linear interpolation:
% Test a scalar function:
seedRNG(6178);
x = (0:10).';  
y = sin(x);
f = chebfun.interp1(x, y, 'linear');
tol = 10*epslevel(f);
pass(1) = norm(feval(f, x) - y) < tol;
pass(2) = numel(f.funs) == 10;
pass(3) = length(f) == 20;

% Test an array-valued function:
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.interp1(x, y, 'linear');
tol = 10*epslevel(f);
pass(4) = norm(feval(f, x) - y) < tol;
pass(5) = numel(f.funs) == 10;
pass(6) = length(f) == 20;

% Test a different domain:
dom = [.01, 10.01];
x = (0:10)';  
y = [sin(x), cos(x)];
f = chebfun.interp1(x, y, 'linear', dom);
err = norm(feval(f, x(2:end-1)) - y(2:end-1,:));
tol = 10*epslevel(f);
pass(7) = all(f.domain == [dom(1), 1:10, dom(2)]);
pass(8) = err < tol;
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

%% Check that things work for row vector input.
x = (0:10).';
y = sin(x);
f1 = chebfun.interp1(x.', y, 'linear');
f2 = chebfun.interp1(x, y.', 'linear');
f3 = chebfun.interp1(x.', y.', 'linear');
pass(21) = (norm(feval(f1, x) - y) < 10*vscale(f1)*epslevel(f1)) && ...
           (norm(feval(f2, x) - y) < 10*vscale(f2)*epslevel(f2)) && ...
           (norm(feval(f3, x) - y) < 10*vscale(f3)*epslevel(f3));
       
% Test random points:
x = rand(11,1);
y = sin(x);
f = chebfun.interp1(x, y, 'linear');
tol = 10*epslevel(f);
pass(22) = norm(feval(f, x) - y) < tol;
pass(23) = numel(f.funs) == 10;
pass(24) = length(f) == 20;

% Test a chebfun
x = chebfun('x', [0, 1]);
y = exp(x);
x = rand(11, 1);
f = chebfun.interp1(x, y, 'poly', [0, 1]);
tol = 100*epslevel(f);
pass(25) = norm(f(x)-y(x), inf) < tol;


%% Test trigonometric interpolation:
x = (0:10)';  
y = cos(pi*x);
f = chebfun.interp1(x, y, 'trig');
tol = 100*epslevel(f);
pass(26) = norm(feval(f, x(2:end-1)) - y(2:end-1)) < tol && ...
    norm(feval(f, x(1))-(y(1)+y(end))/2) < tol ;
pass(27) = numel(f.funs) == 1;
pass(28) = length(f) == length(x)-1;

%%
% Test an array-valued function:
x = (0:10)';  
y = [exp(sin(x)), cos(pi*x)];
f = chebfun.interp1(x, y, 'trig', [0, 11]);
tol = 100*epslevel(f);
pass(29) = norm(feval(f, x) - y) < tol;
pass(30) = numel(f.funs) == 1;
pass(31) = length(f) == length(x);
%%
% Test a different domain:
dom = [-.01, 10.01];
f = chebfun.interp1(x, y, dom);
tol = 100*epslevel(f);
pass(32) = norm(feval(f, x) - y) < tol;
pass(33) = numel(f.funs) == 1;
pass(34) = length(f) == length(x);
pass(35) = all(f.domain == dom);


end

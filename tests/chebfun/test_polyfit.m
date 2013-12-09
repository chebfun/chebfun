% Test file for POLYFIT.

function pass = test_polyfit(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% Test DISCRETEPOLYFIT on [-1 1]:

% Test 1 fits quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
dom = [-1 1];
f = polyfit(chebfun, x, y, 3, dom);
err = norm(feval(f, x) - y);
pass(1) = err < 10*epslevel(f)*vscale(f);

% Test 2 fits cubic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 100)';
y = x.^3;
dom = [-1 1];
f = polyfit(chebfun, x, y, 4, dom);
err = norm(feval(f, x) - y);
pass(2) = err < 5e5*epslevel(f)*vscale(f);

% Test 3 fits quartic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 50)';
y = 2*x + 4*x.^2 + x.^4;
dom = [-1 1];
f = polyfit(chebfun, x, y, 5, dom);
err = norm(feval(f, x) - y);
pass(3) = err < 5e5*epslevel(f)*vscale(f);

% Test 4 fits quartic data {yj} evaluated on Chebyshev points {xj}:
x = chebpts(1000);
y = 2*x + 4*x.^2 + x.^4;
dom = [-1 1];
f = polyfit(chebfun, x, y, 5, dom);
err = norm(feval(f, x) - y);
pass(4) = err < 5e5*epslevel(f)*vscale(f);

% Test 5 fits random data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 10)';
y = -5 + 10*rand(10, 1);
dom = [-1 1];
f = polyfit(chebfun, x, y, 12, dom);
err = norm(feval(f, x) - y);
pass(5) = err < 10*epslevel(f)*vscale(f);

%% Test DISCRETEPOLYFIT on [-1 1] (no domain passed):

% Test 6 fits quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
f = polyfit(chebfun, x, y, 3);
err = norm(feval(f, x) - y);
pass(6) = err < 10*epslevel(f)*vscale(f);

%% Test DISCRETEPOLYFIT on [0 1000] (no domain passed):

% Test 7 fits quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(0, 1000, 1000)';
y = x.^2;
f = polyfit(chebfun, x, y, 3);
err = norm(feval(f, x) - y);
pass(7) = err < 5e1*epslevel(f)*vscale(f);

%% Test CONTINUOUSPOLYFIT on [-1 1]:

% Test 8 checks the 'nothing to do' case with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 6);
err = norm(p-F);
pass(8) = err < epslevel(p)*vscale(p);

% Test 9 fits a quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]);
p = polyfit(F, 5);
err = norm(p-F);
pass(9) = err < epslevel(p)*vscale(p);

% Test 10 fits a quartic chebfun with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 5);
err = norm(p-F);
pass(10) = err < 10*epslevel(p)*vscale(p);

% Test 11 fits a row quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]);
p = polyfit(F.', 5);
err = norm(p-F);
pass(11) = err < epslevel(p)*vscale(p);

%% Test CONTINUOUSPOLYFIT on [0 500]:

% Test 12 fits a quartic chebfun with breakpoint:
F = chebfun(@(x) 3*x.^2, [0 100 500]);
p = polyfit(F, 3);
err = norm(p-F);
pass(12) = err < 1e2*epslevel(p)*vscale(p);

end
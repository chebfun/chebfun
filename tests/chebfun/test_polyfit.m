% Test file for POLYFIT.

function pass = test_polyfit(pref)

if ( nargin == 0 )
    pref = chebpref();
end

seedRNG(7681);

%% Test DISCRETEPOLYFIT on [-1 1]:

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
dom = [-1 1];
f = polyfit(chebfun, x, y, 3, dom);
err = norm(feval(f, x) - y);
pass(1) = err < 10*epslevel(f)*vscale(f);

% Fit cubic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 100)';
y = x.^3;
dom = [-1 1];
f = polyfit(chebfun, x, y, 4, dom);
err = norm(feval(f, x) - y);
pass(2) = err < 5e5*epslevel(f)*vscale(f);

% Fit quartic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 50)';
y = 2*x + 4*x.^2 + x.^4;
dom = [-1 1];
f = polyfit(chebfun, x, y, 5, dom);
err = norm(feval(f, x) - y);
pass(3) = err < 5e5*epslevel(f)*vscale(f);

% Fit quartic data {yj} evaluated on Chebyshev points {xj}:
x = chebpts(1000);
y = 2*x + 4*x.^2 + x.^4;
dom = [-1 1];
f = polyfit(chebfun, x, y, 5, dom);
err = norm(feval(f, x) - y);
pass(4) = err < 5e5*epslevel(f)*vscale(f);

% Fit random data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 10)';
y = -5 + 10*rand(10, 1);
dom = [-1 1];
f = polyfit(chebfun, x, y, 12, dom);
err = norm(feval(f, x) - y);
pass(5) = err < 10*epslevel(f)*vscale(f);

%% Test DISCRETEPOLYFIT on [-1 1] (no domain passed):

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
f = polyfit(chebfun, x, y, 3);
err = norm(feval(f, x) - y);
pass(6) = err < 10*epslevel(f)*vscale(f);

%% Test DISCRETEPOLYFIT on [0 1000] (no domain passed):

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(0, 1000, 1000)';
y = x.^2;
f = polyfit(chebfun, x, y, 3);
err = norm(feval(f, x) - y);
pass(7) = err < 5e1*epslevel(f)*vscale(f);

%% Test DISCRETEPOLYFIT for array-valued inputs:
x = linspace(-1, 1, 10)';
y = [x.^2 x.^3];
f = polyfit(chebfun, x, y, 4);
err = norm(feval(f, x) - y);
pass(8) = err(:) < 10*vscale(f)*epslevel(f);

%% Test CONTINUOUSPOLYFIT on [-1 1]:

% Check the 'nothing to do' case with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 6);
err = norm(p-F);
pass(9) = err < 10*epslevel(p)*vscale(p);

% Fit a quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]);
p = polyfit(F, 5);
err = norm(p-F);
pass(10) = err < 10*epslevel(p)*vscale(p);

% Fit a quartic chebfun with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 5);
err = norm(p-F);
pass(11) = err < 10*epslevel(p)*vscale(p);

% Fit a row quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]).';
p = polyfit(F, 5);
err = norm(p-F);
pass(12) = err < 10*epslevel(p)*vscale(p);

% Fit a row quartic chebfun with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]).';
p = polyfit(F, 5);
err = norm(p-F);
pass(13) = err < 10*epslevel(p)*vscale(p);

%% Test CONTINUOUSPOLYFIT on [0 500]:

% Fit a quartic chebfun with breakpoint:
F = chebfun(@(x) 3*x.^2, [0 100 500]);
p = polyfit(F, 3);
err = norm(p-F);
pass(14) = err < 1e2*epslevel(p)*vscale(p);

%% Test CONTINUOUSPOLYFIT for array-valued inputs:
F = chebfun(@(x) [3*x.^3 (x.^2 + 3*x.^4)], [-1 -0.5 0 0.5 1]);
p = polyfit(F, 5);
err = norm(p-F);
pass(15) = err < 10*epslevel(p)*vscale(p);

end

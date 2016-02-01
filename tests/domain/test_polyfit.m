% Test file for POLYFIT.

function pass = test_polyfit(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(7681);

%% Test POLYFIT on [-1 1]:

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
dom = domain([-1 1]);
f = polyfit(x, y, 2, dom);
err = norm(feval(f, x) - y);
pass(1) = err < 100*eps*vscale(f);

% Fit cubic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 100)';
y = x.^3;
dom = domain([-1 1]);
f = polyfit(x, y, 3, dom);
err = norm(feval(f, x) - y);
pass(2) = err < 5e5*eps*vscale(f);

% Fit quartic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 50)';
y = 2*x + 4*x.^2 + x.^4;
dom = domain([-1 1]);
f = polyfit(x, y, 4, dom);
err = norm(feval(f, x) - y);
pass(3) = err < 5e5*eps*vscale(f);

% Fit quartic data {yj} evaluated on Chebyshev points {xj}:
x = chebpts(1000);
y = 2*x + 4*x.^2 + x.^4;
dom = domain([-1 1]);
f = polyfit(x, y, 4, dom);
err = norm(feval(f, x) - y);
pass(4) = err < 5e5*eps*vscale(f);

% Fit random data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 10)';
y = -5 + 10*rand(10, 1);
dom = domain([-1 1]);
f = polyfit(x, y, 12, dom);
err = norm(feval(f, x) - y);
pass(5) = err < 10*eps*vscale(f);

%% Test POLYFIT on [-1 1]

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(-1, 1, 1000)';
y = x.^2;
dom = domain([-1 1]);
f = polyfit(x, y, 3, dom);
err = norm(feval(f, x) - y);
pass(6) = err < 100*eps*vscale(f);

%% Test POLYFIT on [0 1000] (no domain passed):

% Fit quadratic data {yj} evaluated on equispaced points {xj}:
x = linspace(0, 1000, 1000)';
y = x.^2;
dom = domain([0 1000]);
f = polyfit(x, y, 2, dom);
err = norm(feval(f, x) - y);
pass(7) = err < 5e1*eps*vscale(f);

%% Test POLYFIT for array-valued inputs:
x = linspace(-1, 1, 10)';
y = [x.^2 x.^3];
dom = domain([-1 1]);
f = polyfit(x, y, 3, dom);
err = norm(feval(f, x) - y);
pass(8) = err(:) < 10*vscale(f)*eps;
% Test file for @chebfun/cumsum.m.

function pass = test_cumsum(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty casee.
pass(1) = isempty(cumsum(chebfun()));

% Check an example with breakpoints but no impulses.
f1 = chebfun(@cos, [-1 -0.5 0.5 1], pref);
If1 = cumsum(f1);
If1_exact = @(x) sin(x) - sin(-1);
pass(2) = norm(feval(If1, xr) - If1_exact(xr), inf) < ...
    10*vscale(If1)*epslevel(If1);

% Check behavior for row chebfuns.
f1t = f1.';
If1t = cumsum(f1t);
If1t_exact = @(x) (sin(x) - sin(-1)).';
pass(3) = norm(feval(If1t, xr) - If1t_exact(xr), inf) < ...
    10*vscale(If1t)*epslevel(If1t);

% Check behavior for array-valued chebfuns.
f3 = chebfun(@(x) [cos(x) -sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
If3 = cumsum(f3);
If3_exact = @(x) [sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1);
pass(4) = max(max(abs(feval(If3, xr) - If3_exact(xr)))) < ...
    10*vscale(If3)*epslevel(If3);

f3t = f3.';
If3t = cumsum(f3t);
If3t_exact = @(x) ([sin(x) cos(x) exp(x)] - ...
    repmat([sin(-1) cos(-1) exp(-1)], size(x, 1), 1)).';
pass(5) = max(max(abs(feval(If3t, xr) - If3t_exact(xr)))) < ...
    10*vscale(If3t)*epslevel(If3t);

f4 = chebfun(@(x) [cos(x) -sin(x)], [-1 -0.5 0.5 1], pref);
If4 = cumsum(f4);
If4_exact = @test_If4;
pass(6) = max(max(abs(feval(If4, xr) - If4_exact(xr)))) < ...
    10*vscale(If4)*epslevel(If4);

% Check second argument.
I2f1 = cumsum(f1, 2);
I2f1_exact = @(x) -cos(x) - sin(-1)*x - (-cos(-1) - sin(-1)*-1);
pass(7) = norm(feval(I2f1, xr) - I2f1_exact(xr), inf) < ...
    10*vscale(I2f1)*epslevel(I2f1);

I2f3 = cumsum(f3, 2);
I2f3_exact = @(x) [-cos(x) sin(x) exp(x)] - ...
    [sin(-1)*x cos(-1)*x exp(-1)*x] - ...
    repmat([(-cos(-1) - sin(-1)*-1) (sin(-1) - cos(-1)*-1) ...
            (exp(-1) - exp(-1)*-1)], size(x, 1), 1);
pass(8) = max(max(abs(feval(I2f3, xr) - I2f3_exact(xr)))) < ...
    10*vscale(I2f3)*epslevel(I2f3);

% [TODO]:  Check fractional antiderivatives once implemented.

end

function y = test_If2(x)
    y = zeros(size(x));
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 1);
    y(int1) = sin(x(int1));
    y(int2) = sin(x(int2)) + 1;

    y = y - sin(-1);
end

function y = test_If4(x)
    y = [zeros(size(x)) zeros(size(x))];
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 1);
    y(int1, 1) = sin(x(int1));
    y(int1, 2) = cos(x(int1));
    y(int2, 1) = sin(x(int2));
    y(int2, 2) = cos(x(int2));

    y = y - repmat([sin(-1) cos(-1)], size(x, 1), 1);
end

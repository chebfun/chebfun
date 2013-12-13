% Test file for @chebfun/diff.m.

function pass = test_diff(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty casee.
pass(1) = isempty(diff(chebfun()));

% Check an example with breakpoints:
f1 = chebfun(@sin, [-1 -0.5 0.5 1], pref);
df1 = diff(f1);
df1_exact = @cos;
pass(2) = norm(feval(df1, xr) - df1_exact(xr), inf) < ...
    10*vscale(df1)*epslevel(df1);

% Check behavior for row chebfuns.
f1t = f1.';
df1t = diff(f1t, 1, 2);
df1t_exact = @(x) cos(x).';
pass(3) = norm(feval(df1t, xr) - df1t_exact(xr), inf) < ...
    10*vscale(df1t)*epslevel(df1t);

% Check behavior for array-valued chebfuns.
f3 = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0.5 1], pref);
df3 = diff(f3);
df3_exact = @(x) [cos(x) -sin(x) exp(x)];
pass(4) = max(max(abs(feval(df3, xr) - df3_exact(xr)))) < ...
    10*vscale(df3)*epslevel(df3);

% Check N argument.
d2f1 = diff(f1, 2);
d2f1_exact = @(x) -sin(x);
pass(5) = norm(feval(d2f1, xr) - d2f1_exact(xr), inf) < ...
    1e2*vscale(d2f1)*epslevel(d2f1);

d2f3 = diff(f3, 2);
d2f3_exact = @(x) [-sin(x) -cos(x) exp(x)];
pass(6) = max(max(abs(feval(d2f3, xr) - d2f3_exact(xr)))) < ...
    10*vscale(d2f3)*epslevel(d2f3);

% Check dim argument.
df3_col = diff(f3, 1, 2);
df3_col_exact = @(x) [(cos(x) - sin(x)) (exp(x) - cos(x))];
pass(7) = max(max(abs(feval(df3_col, xr) - df3_col_exact(xr)))) < ...
    10*vscale(df3_col)*epslevel(df3_col);

df3t_row = diff(f3.', 1, 1);
df3t_row_exact = @(x) [(cos(x) - sin(x)) (exp(x) - cos(x))].';
pass(8) = max(max(abs(feval(df3t_row, xr) - df3t_row_exact(xr)))) < ...
    10*vscale(df3t_row)*epslevel(df3t_row);

% Check error conditions.
try
    df3 = diff(f3, 1, 3);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:diff:dim');
end

% [TODO]:  Check fractional derivatives once implemented.

end

function y = test_df2(x)
    y = zeros(size(x));
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 0.5);
    int3 = (0.5 < x) & (x < 1);
    y(int1) = cos(x(int1));
    y(int2) = -sin(x(int2));
    y(int3) = exp(x(int3));
end

function y = test_d2f2(x)
    y = zeros(size(x));
    int1 = (-1 < x) & (x < -0.5);
    int2 = (-0.5 < x) & (x < 0.5);
    int3 = (0.5 < x) & (x < 1);
    y(int1) = -sin(x(int1));
    y(int2) = -cos(x(int2));
    y(int3) = exp(x(int3));
end

function y = test_df4(x)
    y = [zeros(size(x)) zeros(size(x))];
    int1 = (-1 < x(:)) & (x(:) < 0);
    int2 = (0 < x(:)) & (x(:) < 1);
    y(int1, 1) = cos(x(int1));
    y(int1, 2) = -sin(x(int1));
    y(int2, 1) = exp(x(int2));
    y(int2, 2) = exp(x(int2));
end

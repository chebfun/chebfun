% Test file for @chebfun/join.m.

function pass = test_join(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check scalar chebfuns.
f_op = @sin;
f1 = chebfun(f_op, [-1 -0.5 0]);
f2 = chebfun(f_op, [0 0.5 1]);
f = join(f1, f2);
pass(1) = norm(feval(f, xr) - f_op(xr), inf) < 10*vscale(f)*epslevel(f);

% Check array-valued chebfuns.
f_op = @(x) [sin(x) cos(x)];
f1 = chebfun(f_op, [-1 -0.5 0]);
f2 = chebfun(f_op, [0 0.5 1]);
f = join(f1, f2);
err = feval(f, xr) - f_op(xr);
pass(2) = norm(err(:), inf) < 10*vscale(f)*epslevel(f);

% Check row chebfuns.
ft = join(f1.', f2.');
err = feval(ft, xr) - f_op(xr).';
pass(3) = ft.isTransposed && (norm(err(:), inf) < 10*vscale(ft)*epslevel(ft));

% Check operation when the domains don't match.
f = chebfun(@sin, [-1 -0.5 0]);
g = chebfun(@cos, [1 1.5 2]);
h = join(f, g);
pass(4) = isequal(h.domain, [-1 -0.5 0 0.5 1]) && ...
    norm(feval(h, xr) - h_exact(xr), inf) < 10*vscale(h)*epslevel(h);

% Check error conditions.
try
    f = join(f1, f2.')
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:join:trans');
end

% Test on SINGFUN:
op1 = @(x) sin(x)./(x+1);
op2 = @(x) sin(x);
f = chebfun(op1, [-1 -0.5 0], 'exps', [-1 0 0]);
g = chebfun(op2, [1 2]);
h = join(f, g);
x = sort(xr);
h_vals = feval(h, x);
ind = ( x < 0 );
indComp = ~ind;
h_exact1 = feval(op1, x(ind));
h_exact2 = feval(op2, x(indComp)+1);
h_exact = [h_exact1; h_exact2];
err = h_exact - h_vals;
pass(6) = isequal(h.domain, [-1 -0.5 0 1]) && ...
    norm(err, inf) < 5e3*vscale(h)*epslevel(h);

end

function y = h_exact(x)
    y = zeros(size(x));

    xl_ind = x <= 0;
    xl = x(xl_ind);
    y(xl_ind) = sin(xl);

    xr_ind = x > 0;
    xr = x(xr_ind);
    y(xr_ind) = cos(xr + 1);
end

% Test file for @chebfun/addBreaks.m

function pass = test_addBreaks(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Test a scalar chebfun.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = addBreaks(f, [-0.25 0.25]);
pass(1) = isequal(g.domain, [-1 -0.5 -0.25 0 0.25 0.5 1]) && ...
    (norm(feval(f, x) - feval(g, x), Inf) < 10*vscale(f).*epslevel(f));

% Test an array-valued chebfun.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
g = addBreaks(f, [-0.25 0.25]);
err = feval(f, x) - feval(g, x);
pass(2) = isequal(g.domain, [-1 -0.5 -0.25 0 0.25 0.5 1]) && ...
    (norm(err(:), Inf) < 10*vscale(f).*epslevel(f));

end

% Test file for @chebfun/restrict.m

% [TODO]:  Once we're finally handling higher-order impulses, add a test which
% tries restricting a function that has them.

function pass = test_restrict(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Map from [-1, 1] into [a, b] used for transforming test points.
map = @(x, a, b) ((b - a)/2)*(x + 1) + a;

% Test empty inputs.
f = chebfun();
fr = restrict(f, [-0.5 0.5]);
pass(1) = isempty(fr);

f = chebfun(@(x) x, [], pref);
fr = restrict(f, [-0.5 0.5]);
pass(2) = isequal(f.domain, [-1 1]);

% Test behavior when subinterval is same as the original domain.
f_exact = @sin;
f = chebfun(f_exact, [], pref);
fr = restrict(f, [-1 1]);
x = xr;
pass(3) = isequal(fr.domain, [-1 1]) && ...
    (norm(feval(fr, x) - f_exact(x), inf) < 10*fr.vscale*fr.epslevel);

% Test behavior on bogus subinterval input.
try
    fr = restrict(f, [-2 0.5]);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:restrict:subdom');
end

try
    fr = restrict(f, [-0.5 2]);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:restrict:subdom');
end

try
    fr = restrict(f, [-1 -0.25 0.3 0.1 1]);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:restrict:subdom');
end

% Check a scalar function without breakpoints.
f_exact = @(x) 1./(1 + 25*(x - 0.1).^2);
f = chebfun(f_exact, [], pref);

pass(7) = test_restrict_one_function(f, f_exact, [-1 0.5], map, xr);
pass(8) = test_restrict_one_function(f, f_exact, [-0.2 1], map, xr);
pass(9) = test_restrict_one_function(f, f_exact, [-0.2 0.5], map, xr);
pass(10) = test_restrict_one_function(f, f_exact, [-0.35 0.1 0.2 0.5], map, xr);

% Check a scalar function with breakpoints.
f_exact = @(x) 1./(1 + 25*(x - 0.1).^2);
f = chebfun(f_exact, [-1:0.1:1], pref);

pass(11) = test_restrict_one_function(f, f_exact, [-1 0.5], map, xr);
pass(12) = test_restrict_one_function(f, f_exact, [-0.2 1], map, xr);
pass(13) = test_restrict_one_function(f, f_exact, [-0.2 0.5], map, xr);
pass(14) = test_restrict_one_function(f, f_exact, [-0.35 0.1 0.2 0.5], map, xr);

% Check an array-valued function without breakpoints.
f_exact = @(x) [sin(x - 0.1) cos(x + 0.2) exp(x)];
f = chebfun(f_exact, [], pref);

pass(15) = test_restrict_one_function(f, f_exact, [-1 0.5], map, xr);
pass(16) = test_restrict_one_function(f, f_exact, [-0.2 1], map, xr);
pass(17) = test_restrict_one_function(f, f_exact, [-0.2 0.5], map, xr);
pass(18) = test_restrict_one_function(f, f_exact, [-0.35 0.1 0.2 0.5], map, xr);

% Check an array-valued function with breakpoints.
f_exact = @(x) [sin(x - 0.1) cos(x + 0.2) exp(x)];
f = chebfun(f_exact, [-1:0.1:1], pref);

pass(19) = test_restrict_one_function(f, f_exact, [-1 0.5], map, xr);
pass(20) = test_restrict_one_function(f, f_exact, [-0.2 1], map, xr);
pass(21) = test_restrict_one_function(f, f_exact, [-0.2 0.5], map, xr);
pass(22) = test_restrict_one_function(f, f_exact, [-0.35 0.1 0.2 0.5], map, xr);

%% Integration with singfun: piecewise smooth chebfun - splitting on.

% Set the domain:
dom = [-2 7];

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(100*x);
pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
pass(23) = test_restrict_one_function(f, op, dom, map, xr);

end

% Check restriction of a single function.
function pass = test_restrict_one_function(f, f_exact, dom, map, xr)
    fr = restrict(f, dom);
    x = map(xr, dom(1), dom(end));
    err = norm(feval(fr, x) - f_exact(x), inf);
    tol = 10*fr.vscale*fr.epslevel;
    pass = all(ismember(dom, fr.domain)) && ...
        all(err(:) < 2e2*fr.vscale*fr.epslevel);
end

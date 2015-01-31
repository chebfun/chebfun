% Test file for @chebfun/restrict.m

% [TODO]:  Once we're finally handling higher-order impulses, add a test which
% tries restricting a function that has them.

function pass = test_restrict(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
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
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:restrict:subdom');
end

try
    fr = restrict(f, [-0.5 2]);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:restrict:subdom');
end

try
    fr = restrict(f, [-1 -0.25 0.3 0.1 1]);
    pass(6) = false;
catch ME
    pass(6) = true;
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
f = chebfun(f_exact, -1:0.1:1, pref);

pass(19) = test_restrict_one_function(f, f_exact, [-1 0.5], map, xr);
pass(20) = test_restrict_one_function(f, f_exact, [-0.2 1], map, xr);
pass(21) = test_restrict_one_function(f, f_exact, [-0.2 0.5], map, xr);
pass(22) = test_restrict_one_function(f, f_exact, [-0.35 0.1 0.2 0.5], map, xr);

%% Test on singular function: piecewise smooth chebfun - splitting on.

% Set the domain:
dom = [-2 7];
domNew = [-2 1 3.5 6.5 7];
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(100*x).*(x - dom(2)).^pow;
f = chebfun(op, dom, 'exps', [pow pow], 'splitting', 'on');
pass(23) = test_restrict_one_function(f, op, domNew, map, xr);

%% Test on function defined on unbounded domain:

% piecewise function on [-inf b]:

% Set the domain:
dom = [-Inf 1 3*pi];
domRestrict = [-Inf -1 pi 2*pi];
domCheck = [-100 1 2*pi];

% Generate a few random points to use as test values:
x1 = diff(domCheck(1:2)) * rand(100, 1) + domCheck(1);
x2 = diff(domCheck(2:3)) * rand(100, 1) + domCheck(2);

op1 = @(x) exp(x);
op2 = @(x) sin(3*x);
f = chebfun({op1 op2}, dom);
g = restrict(f, domRestrict);
g1Vals = feval(g, x1);
g2Vals = feval(g, x2);
g1Exact = op1(x1);
g2Exact = op2(x2);
err1 = g1Vals - g1Exact;
err2 = g2Vals - g2Exact;
pass(24) = norm([err1 ; err2], inf) < 5*get(g,'epslevel').*get(g,'vscale');

%% Test a bug from issue #528
f = chebfun(@(x) abs(x + 0.04), [-1 0.04 1], 'splitting', 'on');
f = restrict(f, [-0.04 0.04]);
g = chebfun(@(x) abs(x + 0.04), [-0.04 0.04], 'splitting', 'on');
err = norm(f - g, inf);
pass(25) = err < 10*epslevel(g);

%% Test a bug from #727:

f = chebfun(@(x) 4*x.^2-2, [-Inf, Inf]);
g = restrict(f, [-1,1]);
pass(26) = abs(f(1)-g(1)) < epslevel(f)*vscale(f);

%% Test a bug from #1026:

f = chebfun(@(t) t.^0.5./exp(t), [0,inf], 'exps', [0.5 0]);
g = f;
f(1) = f(1); % restrict@unbndfun is called.
pass(27) = ( norm(f-g, inf) < epslevel(f) );

%% Test trigfuns:
f = chebfun(@(x) sin(2*pi*x), [1 2], 'trig');
g = restrict(f, domain(f));
pass(28) = ~isPeriodicTech(g);
g = restrict(f, [1.5, 1.6]);
h = chebfun(@(x) sin(2*pi*x), [1.5, 1.6]);
pass(29) = norm(g-h, inf) < 1e-12;


end

% Check restriction of a single function.
function pass = test_restrict_one_function(f, f_exact, dom, map, xr)
    fr = restrict(f, dom);
    x = map(xr, dom(1), dom(end));
    tol = 10*fr.vscale*fr.epslevel;
    err = norm(feval(fr, x) - f_exact(x), inf);
%     tol2 = 10*fr.vscale*fr.epslevel
%     err2 = max(cellfun(@(d) min(abs(d-fr.domain)), num2cell(dom)))
%     pass = err2 < tol && all(err(:) < tol); % TODO: remove?
    pass = all(ismember(dom, fr.domain)) && ...
        all(err(:) < 5e2*fr.vscale*fr.epslevel); 
end

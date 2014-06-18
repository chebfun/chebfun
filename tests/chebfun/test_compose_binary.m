% Test file for @chebfun/compose.m (binary operators).

function pass = test_compose_binary(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.splitting = 1;

% Smooth operator with smooth functions.
pass(1) = test_one_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) x.*y, pref);

% Smooth operator with smooth functions, non-default domain.
pass(2) = test_one_compose_binary(@(x) cos(2*(x + 0.2)), [-2, 7], ...
    @(x) sin(x - 0.1), [-2 7], @(x, y) x.^2 - y.^2, pref);

% Smooth operator with one non-smooth function input.
pass(3) = test_one_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) sin(x.*y), pref);

% Non-smooth operator with smooth functions, splitting enabled.
pass(4) = test_one_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) abs(x.*y), pref_split);

% Non-smooth operator with one non-smooth function input.
pass(5) = test_one_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) cos(abs(x.*y) - 0.1), pref_split);

% Test array-valued chebfun.
pass(6) = test_one_compose_binary(@(x) [cos(2*(x + 0.2)) sin(2*(x - 0.1))], ...
    [-1 1], @(x) [exp(x) 1./(1 + 25*(x - 0.1).^2)], [-1 1], ...
    @(x, y) x + 2*y, pref);

% Test quasimatrix.
pass(7:9) = test_one_compose_binary_quasi(@(x) [cos(2*(x + 0.2)) sin(2*(x - 0.1))], ...
    [-1 1], @(x) [exp(x) 1./(1 + 25*(x - 0.1).^2)], [-1 1], ...
    @(x, y) x + 2*y, pref);

%% Test for singular function:
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

f = chebfun(@(x) sin(20*x)./(x-dom(1)), dom, 'exps', [-1 0]);
g = chebfun(@(x) x.^2./(x-dom(1)), dom, 'exps', [-1 0]);
h = compose(f, @plus, g);
hVals = feval(h, x);
oph = @(x) (sin(20*x)+x.^2)./(x-dom(1));
hExact = oph(x);
err = hVals - hExact;
pass(10) = ( norm(err, inf) < 1e5*vscale(h)*epslevel(h) ) && ...
    ( isequal(oph(dom(2)), feval(h, dom(2))) );


%% Compose two functions defined on an unbounded domain, i.e. F + G:

% Set the domain:
dom = [0 Inf];
domCheck = [0 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x);
opg = @(x) x.*exp(-x);
oph = @(x) (x+1).*exp(-x);
f = chebfun(opf, dom);
g = chebfun(opg, dom);
h = compose(f, @plus, g);
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(11) = norm(err, inf) < get(h,'epslevel')*get(h,'vscale');

end

% Test composition with binary operators.
function pass = test_one_compose_binary(f_exact, fdom, g_exact, gdom, op, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = chebfun(f_exact, fdom, pref);
    g = chebfun(g_exact, gdom, pref);
    h = compose(f, op, g, pref);
    h_exact = @(x) op(f_exact(x), g_exact(x));
    dom = union(fdom, gdom);
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(h, x) - h_exact(x), inf);
    pass = (err < 20*vscale(h)*epslevel(h)) && ...
        isequal(h_exact(dom), feval(h, dom));
end

% Test composition with binary operators.
function pass = test_one_compose_binary_quasi(f_exact, fdom, g_exact, gdom, op, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = chebfun(f_exact, fdom, pref);
    g = chebfun(g_exact, gdom, pref);
    fq = quasimatrix(f_exact, fdom, pref);
    gq = quasimatrix(g_exact, gdom, pref); 
    
    h_exact = @(x) op(f_exact(x), g_exact(x));
    dom = union(fdom, gdom);
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    
    h = compose(f, op, gq, pref);
    err = norm(feval(h, x) - h_exact(x), inf);
    pass(1) = (err < 20*vscale(h)*epslevel(h)) && ...
        isequal(h_exact(dom), feval(h, dom));
    
    h = compose(fq, op, g, pref);
    err = norm(feval(h, x) - h_exact(x), inf);
    pass(2) = (err < 20*vscale(h)*epslevel(h)) && ...
        isequal(h_exact(dom), feval(h, dom));
    
    h = compose(fq, op, gq, pref);
    err = norm(feval(h, x) - h_exact(x), inf);
    pass(3) = (err < 20*vscale(h)*epslevel(h)) && ...
        isequal(h_exact(dom), feval(h, dom));
end

% Test file for @chebfun/compose.m (binary operators).

function pass = test_compose(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.enableBreakpointDetection = 1;

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
    pass = (err < 20*h.vscale*h.epslevel) && ...
        isequal(h_exact(dom), feval(h, dom));
end

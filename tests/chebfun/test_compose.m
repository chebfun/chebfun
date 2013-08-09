% Test file for @chebfun/compose.m

function pass = test_compose(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.chebfun.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Test empty input.
f = chebfun();
g = compose(f, @sin, [], pref);
pass(1) = isempty(g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test composition with unary operators.

% Smooth operator with smooth function.
pass(2) = test_compose_unary(@(x) cos(x - 0.1), [-1 1], @exp, pref);

% Smooth operator with smooth function, non-default domain.
pass(3) = test_compose_unary(@(x) 1./(1 + (x - 0.1).^2), [-5 5], ...
    @(x) sin(20*x), pref);

% Smooth operator with non-smooth function.
pass(4) = test_compose_unary(@(x) abs(x + 0.2) + abs(x - 0.3), ...
    [-1 -0.2 0.3 1], @(x) cos(x.^2), pref);

% Non-smooth operator with smooth function, splitting disabled.
warnstate = warning('off');
test_compose_unary(@(x) sin(10*(x - 0.1)), [-1 1], @abs, pref);
[warnmsg, warnid] = lastwarn();
pass(5) = strcmp(warnid, 'CHEBFUN:compose:resolve');
warning(warnstate);

% Non-smooth operator with smooth function, splitting enabled.
pass(6) = test_compose_unary(@(x) sin(10*(x - 0.1)), [-1 1], @abs, pref_split);

% Non-smooth operator with non-smooth function.
pass(7) = test_compose_unary(@(x) abs(sin(3*(x - 0.1))) - 0.5, ...
    [-1 -0.95 0.1 1], @abs, pref_split);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test composition with binary operators.

% Smooth operator with smooth functions.
pass(8) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) x.*y, pref);

% Smooth operator with smooth functions, non-default domain.
pass(9) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-2, 7], ...
    @(x) sin(x - 0.1), [-2 7], @(x, y) x.^2 - y.^2, pref);

% Smooth operator with one non-smooth function input.
pass(10) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) sin(x.*y), pref);

% Non-smooth operator with smooth functions, splitting enabled.
pass(11) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) abs(x.*y), pref_split);

% Non-smooth operator with one non-smooth function input.
pass(12) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) cos(abs(x.*y) - 0.1), pref_split);

% [TODO]:  Once composeChebfuns() is working, add tests for it.

end

% Test composition of a function with a unary operator.
function pass = test_compose_unary(f_exact, dom, op, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = chebfun(f_exact, dom, pref);
    g = compose(f, op, [], pref);
    g_exact = @(x) op(f_exact(x));
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(g, x) - g_exact(x), inf);
    pass = (err < 20*g.vscale*g.epslevel) && ...
        isequal(g_exact(f.domain), feval(g, f.domain));
end

% Test composition with binary operators.
function pass = test_compose_binary(f_exact, fdom, g_exact, gdom, op, pref)
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

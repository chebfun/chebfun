% Test file for @chebfun/compose.m (unary operators.

function pass = test_compose_unary(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.chebfun.splitting = 1;

% Test empty input.
f = chebfun();
g = compose(f, @sin, [], pref);
pass(1) = isempty(g);

% Smooth operator with smooth function.
pass(2) = test_one_compose_unary(@(x) cos(x - 0.1), [-1 1], @exp, pref);

% Smooth operator with smooth function, non-default domain.
pass(3) = test_one_compose_unary(@(x) 1./(1 + (x - 0.1).^2), [-5 5], ...
    @(x) sin(20*x), pref);

% Smooth operator with non-smooth function.
pass(4) = test_one_compose_unary(@(x) abs(x + 0.2) + abs(x - 0.3), ...
    [-1 -0.2 0.3 1], @(x) cos(x.^2), pref);

% Non-smooth operator with smooth function, splitting disabled.
warnstate = warning('off');
test_one_compose_unary(@(x) sin(10*(x - 0.1)), [-1 1], @abs, pref);
[warnmsg, warnid] = lastwarn();
pass(5) = strcmp(warnid, 'CHEBFUN:compose:resolve');
warning(warnstate);

% Non-smooth operator with smooth function, splitting enabled.
pass(6) = test_one_compose_unary(@(x) sin(10*(x - 0.1)), [-1 1], @abs, ...
    pref_split);

% Non-smooth operator with non-smooth function.
pass(7) = test_one_compose_unary(@(x) abs(sin(3*(x - 0.1))) - 0.5, ...
    [-1 -0.95 0.1 1], @abs, pref_split);

% Test array-valued chebfun.
pass(8) = test_one_compose_unary(@(x) [cos(2*(x + 0.2)), sin(2*(x - 0.1))], ...
    [-1 1], @exp, pref);

end

% Test composition of a function with a unary operator.
function pass = test_one_compose_unary(f_exact, dom, op, pref)
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

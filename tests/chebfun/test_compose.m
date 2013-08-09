% Test file for @chebfun/compose.m

function pass = test_compose(pref)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Test array-valued chebfun.
pass(8) = test_compose_unary(@(x) [cos(2*(x + 0.2)), sin(2*(x - 0.1))], ...
    [-1 1], @exp, pref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test composition with binary operators.

% Smooth operator with smooth functions.
pass(9) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) x.*y, pref);

% Smooth operator with smooth functions, non-default domain.
pass(10) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-2, 7], ...
    @(x) sin(x - 0.1), [-2 7], @(x, y) x.^2 - y.^2, pref);

% Smooth operator with one non-smooth function input.
pass(11) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) sin(x.*y), pref);

% Non-smooth operator with smooth functions, splitting enabled.
pass(12) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) sin(x - 0.1), [-1 1], @(x, y) abs(x.*y), pref_split);

% Non-smooth operator with one non-smooth function input.
pass(13) = test_compose_binary(@(x) cos(2*(x + 0.2)), [-1, 1], ...
    @(x) abs(x - 0.1), [-1 0.1 1], @(x, y) cos(abs(x.*y) - 0.1), pref_split);

% Test array-valued chebfun.
pass(14) = test_compose_binary(@(x) [cos(2*(x + 0.2)) sin(2*(x - 0.1))], ...
    [-1 1], @(x) [exp(x) 1./(1 + 25*(x - 0.1).^2)], [-1 1], ...
    @(x, y) x + 2*y, pref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test composition of two chebfuns.

% [TODO]:  Once composeChebfuns() can handle domain checks, test that too.

% Two smooth chebfuns.
pass(15) = test_compose_chebfuns(@(x) cos(2*(x + 0.2)), @(x) sin(x - 0.1), ...
    [-1 1], pref);

% Two smooth chebfuns, non-default domain.
pass(16) = test_compose_chebfuns(@(x) cos(2*(x + 0.2)), @(x) sin(x - 0.1), ...
    [-2 7], pref);

% Smooth chebfun of a non-smooth chebfun.
pass(17) = test_compose_chebfuns(@(x) abs(cos(2*(x + 0.2))), ...
    @(x) sin(x - 0.1), [-1 1], pref_split);

% Non-smooth chebfun of a smooth chebfun.
pass(18) = test_compose_chebfuns(@(x) sin(x - 0.1), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Non-smooth chebfun of a non-smooth chebfun.
pass(19) = test_compose_chebfuns(@(x) abs(sin(2*(x - 0.1))), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Can't compose two array-valued chebfuns.
try
    test_compose_chebfuns(@(x) [sin(x) cos(x)], @(x) [exp(x) -sin(x)], ...
        [-1 1], pref);
    pass(20) = false;
catch ME
    pass(20) = true;
end

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

% Test composition of two chebfuns.
function pass = test_compose_chebfuns(f_exact, g_exact, dom, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = chebfun(f_exact, dom, pref);
    g = chebfun(g_exact, dom, pref);
    h = compose(f, g, [], pref);
    h_exact = @(x) g_exact(f_exact(x));
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(h, x) - h_exact(x), inf);
    pass = (err < 20*h.vscale*h.epslevel) && ...
        isequal(feval(g, feval(f, (dom))), feval(h, dom));
end

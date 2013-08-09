% Test file for @chebfun/compose.m (composition of chebfuns).

function pass = test_compose(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.chebfun.splitting = 1;

% [TODO]:  Once composeChebfuns() can handle domain checks, test that too.

% Two smooth chebfuns.
pass(1) = test_one_compose_chebfuns(@(x) cos(2*(x + 0.2)), ...
    @(x) sin(x - 0.1), [-1 1], pref);

% Two smooth chebfuns, non-default domain.
pass(2) = test_one_compose_chebfuns(@(x) cos(2*(x + 0.2)), ...
    @(x) sin(x - 0.1), [-2 7], pref);

% Smooth chebfun of a non-smooth chebfun.
pass(3) = test_one_compose_chebfuns(@(x) abs(cos(2*(x + 0.2))), ...
    @(x) sin(x - 0.1), [-1 1], pref_split);

% Non-smooth chebfun of a smooth chebfun.
pass(4) = test_one_compose_chebfuns(@(x) sin(x - 0.1), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Non-smooth chebfun of a non-smooth chebfun.
pass(5) = test_one_compose_chebfuns(@(x) abs(sin(2*(x - 0.1))), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Can't compose two array-valued chebfuns.
try
    test_one_compose_chebfuns(@(x) [sin(x) cos(x)], @(x) [exp(x) -sin(x)], ...
        [-1 1], pref);
    pass(6) = false;
catch ME
    pass(6) = true;
end

end

% Test composition of two chebfuns.
function pass = test_one_compose_chebfuns(f_exact, g_exact, dom, pref)
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


% Test file for @chebfun/compose.m (unary operators.

function pass = test_compose_unary(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.splitting = 1;

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
pass(5) = strcmp(warnid, 'CHEBFUN:CHEBFUN:compose:resolve');
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

% Test quasimatrix.
pass(9) = test_one_compose_unary_quasi(@(x) [cos(2*(x + 0.2)), ...
    sin(2*(x - 0.1))], [-1 1], @exp, pref);

%% Test for singular function:
dom = [0 1];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

f = chebfun(@(x) sqrt(x), dom, 'blowup', 2);
g = sin(f);
gVals = feval(g, x);
opg = @(x) sin(sqrt(x));
gExact = opg(x);
err = gVals - gExact;
pass(10) = ( norm(err, inf) < 1e1*vscale(g)*eps ) && ...
        isequal(opg(domain(f)), feval(g, domain(f)));
    
%% Compose a function defined on an unbounded domain with an operator, i.e. 
% OP(F)

% Set the domain:
dom = [0 Inf];
domCheck = [0 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x);
opg = @(x) sin(exp(-x));
f = chebfun(opf, dom);
g = compose(f, @sin);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(11) = norm(err, inf) < 1e1*eps*get(g,'vscale');


%% Test related to #686

vals = [...
   3.141592652798157
   3.141592652796547
   3.141592652772516
   3.141592652667793
   3.141592652373061
   3.141592651675060
   3.141592650102323
   3.141592646526151
   3.141592638146115
   3.141592617916730
   3.141592568136906
   3.141592445124339
   3.141592145023115
   3.141591435222665
   3.141589837802166
   3.141586481864910
   3.141580027749437
   3.141568894324190
   3.141552052732393
   3.141530333262496
   3.141507475435311
   3.141489646225016
   3.141482861866948];
u = chebfun(vals, 'tech', @chebtech2);
v = sin(u);
pass(12) = ishappy(v) && length(v) < 1e3;

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
    pass = (err < 1e2*vscale(g)*eps) && ...
        isequal(g_exact(domain(f)), feval(g, domain(f)));
    
end

% Test composition of a function with a unary operator.
function pass = test_one_compose_unary_quasi(f_exact, dom, op, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = quasimatrix(f_exact, dom, pref);
    g = compose(f, op, [], pref);
    g_exact = @(x) op(f_exact(x));
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(g, x) - g_exact(x), inf);
    pass = (err < 20*vscale(g)*eps) && ...
        isequal(g_exact(domain(f)), feval(g, domain(f)));
end

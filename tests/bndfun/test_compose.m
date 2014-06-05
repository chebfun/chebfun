% Test file for bndfun/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set the domain
dom = [-2 7];
x = linspace(dom(1), dom(2), 1000);

% Compose a scalar-valued BNDFUN object with sin(x):
f = bndfun(@(x) x, struct('domain', dom));
g = compose(f, @sin, [], [], pref);
h = @(x) sin(x);
pass(1) = norm(h(x) - feval(g, x), inf) < ...
    10*max(get(g, 'vscale'))*get(g, 'epslevel');

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x], struct('domain', dom));
g = compose(f, @sin, [], [], pref);
h = @(x) [sin(x) sin(x)];
err = feval(g, x) - h(x);
pass(2) = norm(err(:), inf) < 10*max(get(g, 'vscale').*get(g, 'epslevel'));

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x.^2], struct('domain', dom));
g = compose(f, @sin, [], [], pref);
pass(3) = norm(sin([x, x.^2]) - feval(g, x), inf) < ...
    10*max(get(g,'vscale').*get(g, 'epslevel'));

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x x.^2], struct('domain', dom));
g = compose(f, @sin, [], [], pref);
pass(4) = norm(sin([x x x.^2]) - feval(g, x), inf) < ...
    10*max(get(g,'vscale').*get(g, 'epslevel'));

% Compose 2 BNDFUN objects with a binary function:
f1 = bndfun(@(x) sin(x), struct('domain', dom));
f2 = bndfun(@(x) cos(x), struct('domain', dom));
g = compose(f1, @plus, f2, [], pref);
h = @(x) sin(x) + cos(x);
pass(5) = norm(h(x) - feval(g, x), inf) < ...
    10*max(get(g, 'vscale').*get(g, 'epslevel'));

% Compose 2 array-valued BNDFUN objects with a binary function:
f1 = bndfun(@(x) [sin(x) cos(x)], struct('domain', dom));
f2 = bndfun(@(x) [cos(x) exp(x)], struct('domain', dom));
g = compose(f1, @times, f2, [], pref);
h = bndfun(@(x) [sin(x).*cos(x) cos(x).*exp(x)], struct('domain', dom));
pass(6) = norm(feval(h, x) - feval(g, x)) < ...
    10*max(get(g,'vscale').*get(g, 'epslevel'));

% Compose f(g), when f and g are BNDFUN objects:
f = bndfun(@(x) x.^2, struct('domain', dom));
g = bndfun(@(x) sin(x), struct('domain', [0 dom(2)^2]));
h = compose(f, g);
pass(7) = norm(feval(h, x) - sin(x.^2), inf) < ...
    10*max(get(h,'vscale').*get(g, 'epslevel'));

% Compose f(g), when f and g are BNDFUN objects and g is array-valued:
f = bndfun(@(x) x.^2, struct('domain', dom));
g = bndfun(@(x) [sin(x) cos(x)], struct('domain', [0 dom(2)^2]));
h = compose(f, g);
pass(8) = norm(feval(h, x) - [sin(x.^2) cos(x.^2)], inf) < ...
    10*max(get(h,'vscale').*get(h, 'epslevel'));

% Compose f(g), when f and g are BNDFUN objects and f is array-valued:
f = bndfun(@(x) [x x.^2], struct('domain', dom));
g = bndfun(@(x) sin(x), struct('domain', [dom(1) dom(2)^2]));
h = compose(f, g);
pass(9) = norm(feval(h, x) - [sin(x) sin(x.^2)], inf) < ...
    10*max(get(h,'vscale').*get(h, 'epslevel'));

end

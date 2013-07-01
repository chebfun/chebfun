% Test file for bndfun/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = fun.pref;
    pref = chebtech.pref(pref);
end

% Set the tolerance
tol = 10*pref.fun.eps;

% Set the domain
dom = [-2 7];
x = linspace(dom(1), dom(2), 1000);

pass = zeros(1, 12); % Pre-allocate pass matrix.

% Compose a scalar-valued BNDFUN object with sin(x):
f = bndfun(@(x) x, dom);
g = compose(f, @sin, pref);
h = bndfun(@sin, dom);
pass(1) = normest(h - g) < tol;

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x], dom);
g = compose(f, @sin, pref);
h = bndfun(@(x) [sin(x) sin(x)], dom);
pass(2) = normest(h - g) < tol;

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x.^2], dom);
g = compose(f, @sin, pref);

pass(3) = norm(sin([x, x.^2]) - feval(g, x), inf) < ...
    max(g.onefun.vscale)*tol;

% Compose an array-valued BNDFUN object with sin(x):
f = bndfun(@(x) [x x x.^2], dom);
g = compose(f, @sin, pref);
pass(4) = norm(sin([x x x.^2]) - feval(g, x), inf) < ...
    max(g.onefun.vscale)*tol;

% Compose 2 BNDFUN objects with a binary function:
f1 = bndfun(@(x) sin(x), dom);
f2 = bndfun(@(x) cos(x), dom);
g = compose(f1, @plus, f2, pref);
h = @(x) sin(x) + cos(x);
pass(5) = norm(h(x) - feval(g, x), inf) < 5*tol;

% Compose 2 array-valued BNDFUN objects with a binary function:
f1 = bndfun(@(x) [sin(x) cos(x)], dom);
f2 = bndfun(@(x) [cos(x) exp(x)], dom);
g = compose(f1, @times, f2, pref);
h = bndfun(@(x) [sin(x).*cos(x) cos(x).*exp(x)], dom);
pass(6) = normest(h - g) < max(g.onefun.vscale)*tol;

% Compose f(g), when f and g are BNDFUN objects:
f = bndfun(@(x) x.^2, dom);
g = bndfun(@(x) sin(x), [0 dom(2)^2]);
h = compose(f, g);
pass(7) = norm(feval(h, x) - sin(x.^2), inf) < max(h.onefun.vscale)*tol;

% Compose f(g), when f and g are BNDFUN objects and g is array-valued:
f = bndfun(@(x) x.^2, dom);
g = bndfun(@(x) [sin(x) cos(x)], [0 dom(2)^2]);
h = compose(f, g);
pass(8) = norm(feval(h, x) - [sin(x.^2) cos(x.^2)], inf) < ...
    max(h.onefun.vscale)*tol;

% Compose f(g), when f and g are BNDFUN objects and f is array-valued:
f = bndfun(@(x) [x x.^2], dom);
g = bndfun(@(x) sin(x), [dom(1) dom(2)^2]);
h = compose(f, g);
pass(9) = norm(feval(h, x) - [sin(x) sin(x.^2)], inf) < ...
    max(h.onefun.vscale)*tol;

% We cannot expect to compose two array-valued BNDFUN objects f(g):
try
    f = bndfun(@(x) [x x.^2], dom);
    g = bndfun(@(x) [sin(x), cos(x)], dom);
    compose(f, g);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'BNDFUN:compose:domainMismatch');
end

% We cannot expect to compose two array-valued BNDFUN objects in this way:
try
    f = bndfun(@(x) [x x.^2], dom);
    g = bndfun(@(x) sin(x), dom);
    compose(f, @plus, g)
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:dim');
end

% Can't compose two BNDFUN objects f(g) if the range of g does not lie in
% the domain of f:
try
    f = bndfun(@(x) sin(x), dom);
    g = bndfun(@(x) 100*cos(x), dom);
    compose(g, f);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'BNDFUN:compose:domainMismatch');
end

end


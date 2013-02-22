function pass = test_compose(pref)

if ( nargin < 1 )
    pref = funcheb2.pref;
end
tol = 10*pref.funcheb2.eps;

% % Compose a scalar-valued FUNCHEB2 object with sin(x):
f = funcheb2(@(x) x);
g = compose(f, @sin, pref);
h = funcheb2(@sin);
pass(1) = norm(h.values - g.values, inf) < tol;

% Compose a multi-valued FUNCHEB2 object with sin(x):
f = funcheb2(@(x) [x x]);
g = compose(f, @sin, pref);
h = funcheb2(@(x) [sin(x) sin(x)]);
pass(2) = norm(h.values - g.values, inf) < tol;

% Compose a multi-valued FUNCHEB2 object with sin(x):
f = funcheb2(@(x) [x x.^2]);
g = compose(f, @sin, pref);
x = funcheb2.chebpts(length(g));
pass(3) = norm(sin([x, x.^2])- g.values, inf) < tol;

% Compose a multivalued FUNCHEB2 object with sin(x):
f = funcheb2(@(x) [x x x.^2]);
g = compose(f, @sin, pref);
x = funcheb2.chebpts(length(g));
pass(4) = norm(sin([x x x.^2]) - g.values, inf) < tol;

% Compose 2 FUNCHEB2 objects with a binary function:
f1 = funcheb2(@(x) sin(x));
f2 = funcheb2(@(x) cos(x));
g = compose(f1, @plus, f2, pref);
h = funcheb2(@(x) sin(x) + cos(x));
pass(5) = norm(h.values - g.values, inf) < tol;

% Compose 2 multivalued FUNCHEB2 objects with a binary function:
f1 = funcheb2(@(x) [sin(x) cos(x)]);
f2 = funcheb2(@(x) [cos(x) exp(x)]);
g = compose(f1, @times, f2, pref);
h = funcheb2(@(x) [sin(x).*cos(x) cos(x).*exp(x)]);
pass(6) = norm(h.values - g.values, inf) < tol;

% Compose f(g), when f and g are FUNCHEB2 objects:
f = funcheb2(@(x) x.^2);
g = funcheb2(@(x) sin(x));
h = compose(f, g);
x = funcheb2.chebpts(length(h));
pass(7) = norm(h.values - sin(x.^2), inf) < tol;

% Compose f(g), when f and g are FUNCHEB2 objects and g is multivalued:
f = funcheb2(@(x) x.^2);
g = funcheb2(@(x) [sin(x) cos(x)]);
h = compose(f, g);
x = funcheb2.chebpts(length(h));
pass(8) = norm(h.values - [sin(x.^2) cos(x.^2)], inf) < tol;

% Compose f(g), when f and g are FUNCHEB2 objects and f is multivalued:
f = funcheb2(@(x) [x x.^2]);
g = funcheb2(@(x) sin(x));
h = compose(f, g);
x = funcheb2.chebpts(length(h));
pass(9) = norm(h.values - [sin(x) sin(x.^2)], inf) < tol;

% We cannot expect to compose two multivalued FUNCHEB2 objects f(g):
try 
    f = funcheb2(@(x) [x x.^2]);
    g = funcheb2(@(x) [sin(x), cos(x)]);
    compose(f, g);
    pass(10) = false;
catch ME 
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:multival');
end

% We cannot expect to compose two multivalued FUNCHEB2 objects in this way:
try 
    f = funcheb2(@(x) [x x.^2]);
    g = funcheb2(@(x) sin(x));
    compose(f, @plus, g)
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:dim');
end

try
    f = funcheb2(@(x) sin(x));
    g = funcheb2(@(x) 100*cos(x));
    compose(g,f);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:range');
end

end

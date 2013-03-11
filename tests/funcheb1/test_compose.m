function pass = test_compose(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

% % Compose a scalar-valued FUNCHEB1 object with sin(x):
f = funcheb1(@(x) x);
g = compose(f, @sin, pref);
h = funcheb1(@sin);
pass(1) = norm(h.values - g.values, inf) < tol;

% Compose a multi-valued FUNCHEB1 object with sin(x):
f = funcheb1(@(x) [x x]);
g = compose(f, @sin, pref);
h = funcheb1(@(x) [sin(x) sin(x)]);
pass(2) = norm(h.values - g.values, inf) < tol;

% Compose a multi-valued FUNCHEB1 object with sin(x):
f = funcheb1(@(x) [x x.^2]);
g = compose(f, @sin, pref);
x = funcheb1.chebpts(length(g));
pass(3) = norm(sin([x, x.^2])- g.values, inf) < tol;

% Compose a multi-valued FUNCHEB1 object with sin(x):
f = funcheb1(@(x) [x x x.^2]);
g = compose(f, @sin, pref);
x = funcheb1.chebpts(length(g));
pass(4) = norm(sin([x x x.^2]) - g.values, inf) < tol;

% Compose 2 FUNCHEB1 objects with a binary function:
f1 = funcheb1(@(x) sin(x));
f2 = funcheb1(@(x) cos(x));
g = compose(f1, @plus, f2, pref);
h = funcheb1(@(x) sin(x) + cos(x));
pass(5) = norm(h.values - g.values, inf) < tol;

% Compose 2 multivalued FUNCHEB1 objects with a binary function:
f1 = funcheb1(@(x) [sin(x) cos(x)]);
f2 = funcheb1(@(x) [cos(x) exp(x)]);
g = compose(f1, @times, f2, pref);
h = funcheb1(@(x) [sin(x).*cos(x) cos(x).*exp(x)]);
pass(6) = norm(h.values - g.values, inf) < tol;

% Compose g(f), when f and g are FUNCHEB1 objects:
f = funcheb1(@(x) x.^2);
g = funcheb1(@(x) sin(x));
h = compose(f, g);
x = funcheb1.chebpts(length(h));
pass(7) = norm(h.values - sin(x.^2), inf) < tol;

% Compose g(f), when f and g are FUNCHEB1 objects and g is multivalued:
f = funcheb1(@(x) x.^2);
g = funcheb1(@(x) [sin(x) cos(x)]);
h = compose(f, g);
x = funcheb1.chebpts(length(h));
pass(8) = norm(h.values - [sin(x.^2) cos(x.^2)], inf) < tol;

% Compose g(f), when f and g are FUNCHEB1 objects and f is multivalued:
f = funcheb1(@(x) [x x.^2]);
g = funcheb1(@(x) sin(x));
h = compose(f, g);
x = funcheb1.chebpts(length(h));
pass(9) = norm(h.values - [sin(x) sin(x.^2)], inf) < tol;

% We cannot expect to compose two multivalued FUNCHEB1 objects g(f):
try 
    f = funcheb1(@(x) [x x.^2]);
    g = funcheb1(@(x) [sin(x), cos(x)]);
    compose(f, g);
    pass(10) = false;
catch ME 
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:multival');
end

% We cannot expect to compose two multivalued FUNCHEB1 objects with a 
% binary function:
try 
    f = funcheb1(@(x) [x x.^2]);
    g = funcheb1(@(x) sin(x));
    compose(f, @plus, g)
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:dim');
end

% The range of f must coincide with the domain of g for the composition g(f)
try
    f = funcheb1(@(x) 100*cos(x));
    g = funcheb1(@(x) sin(x));
    compose(f,g);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:range');
end

end

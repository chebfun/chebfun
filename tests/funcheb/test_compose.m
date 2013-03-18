% Test file for funcheb/compose.

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    tol = 10*pref.funcheb.eps;

    % Compose a scalar-valued FUNCHEB object with sin(x):
    f = testclass.make(@(x) x);
    g = compose(f, @sin, pref);
    h = testclass.make(@sin);
    pass(n, 1) = norm(h.values - g.values, inf) < tol;

    % Compose a multi-valued FUNCHEB object with sin(x):
    f = testclass.make(@(x) [x x]);
    g = compose(f, @sin, pref);
    h = testclass.make(@(x) [sin(x) sin(x)]);
    pass(n, 2) = norm(h.values - g.values, inf) < tol;

    % Compose a multi-valued FUNCHEB object with sin(x):
    f = testclass.make(@(x) [x x.^2]);
    g = compose(f, @sin, pref);
    x = g.points();
    pass(n, 3) = norm(sin([x, x.^2])- g.values, inf) < tol;
    
    % Compose a multivalued FUNCHEB object with sin(x):
    f = testclass.make(@(x) [x x x.^2]);
    g = compose(f, @sin, pref);
    x = g.points();
    pass(n, 4) = norm(sin([x x x.^2]) - g.values, inf) < tol;
    
    % Compose 2 FUNCHEB objects with a binary function:
    f1 = testclass.make(@(x) sin(x));
    f2 = testclass.make(@(x) cos(x));
    g = compose(f1, @plus, f2, pref);
    h = testclass.make(@(x) sin(x) + cos(x));
    pass(n, 5) = norm(h.values - g.values, inf) < tol;
    
    % Compose 2 multivalued FUNCHEB objects with a binary function:
    f1 = testclass.make(@(x) [sin(x) cos(x)]);
    f2 = testclass.make(@(x) [cos(x) exp(x)]);
    g = compose(f1, @times, f2, pref);
    h = testclass.make(@(x) [sin(x).*cos(x) cos(x).*exp(x)]);
    pass(n, 6) = norm(h.values - g.values, inf) < tol;
    
    % Compose f(g), when f and g are FUNCHEB objects:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 7) = norm(h.values - sin(x.^2), inf) < tol;
    
    % Compose f(g), when f and g are FUNCHEB objects and g is multivalued:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) [sin(x) cos(x)]);
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 8) = norm(h.values - [sin(x.^2) cos(x.^2)], inf) < tol;
    
    % Compose f(g), when f and g are FUNCHEB objects and f is multivalued:
    f = testclass.make(@(x) [x x.^2]);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 9) = norm(h.values - [sin(x) sin(x.^2)], inf) < tol;
    
    % We cannot expect to compose two multivalued FUNCHEB objects f(g):
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) [sin(x), cos(x)]);
        compose(f, g);
        pass(n, 10) = false;
    catch ME 
        pass(n, 10) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:multival');
    end
    
    % We cannot expect to compose two multivalued FUNCHEB objects in this way:
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) sin(x));
        compose(f, @plus, g)
        pass(n, 11) = false;
    catch ME
        pass(n, 11) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:dim');
    end
    
    try
        f = testclass.make(@(x) sin(x));
        g = testclass.make(@(x) 100*cos(x));
        compose(g,f);
        pass(n, 12) = false;
    catch ME
        pass(n, 12) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:compose:range');
    end
end

end

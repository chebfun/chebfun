% Test file for chebtech/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:4
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();

        if ( n == 3 )
            pref.refinementFunction = 'resampling';
        end

        if ( n == 4 )
            pref.extrapolate = true;
        end
    end

    % Compose a scalar-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) x);
    g = compose(f, @sin, [], pref);
    h = testclass.make(@sin);
    pass(n, 1) = norm(h.values - g.values, inf) < 10*h.vscale.*h.epslevel;

    % Compose an array-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x]);
    g = compose(f, @sin, [], pref);
    h = testclass.make(@(x) [sin(x) sin(x)]);
    pass(n, 2) = norm(h.values - g.values, inf) < ...
        10*max(h.vscale.*h.epslevel);

    % Compose an array-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x.^2]);
    g = compose(f, @sin, [], pref);
    x = g.points();
    pass(n, 3) = norm(sin([x, x.^2])- g.values, inf) < ...
        10*max(h.vscale.*h.epslevel);
    
    % Compose an array-valued CHEBTECH object with sin(x):
    f = testclass.make(@(x) [x x x.^2]);
    g = compose(f, @sin, [], pref);
    x = g.points();
    pass(n, 4) = norm(sin([x x x.^2]) - g.values, inf) < ...
        10*max(h.vscale.*h.epslevel);
    
    % Compose 2 CHEBTECH objects with a binary function:
    f1 = testclass.make(@(x) sin(x));
    f2 = testclass.make(@(x) cos(x));
    g = compose(f1, @plus, f2, pref);
    x = g.points;
    h = testclass.make(sin(x) + cos(x));
    pass(n, 5) = norm(h.values - g.values, inf) < 10*h.vscale.*h.epslevel;
    
    % Compose 2 array-valued CHEBTECH objects with a binary function:
    f1 = testclass.make(@(x) [sin(x) cos(x)]);
    f2 = testclass.make(@(x) [cos(x) exp(x)]);
    g = compose(f1, @times, f2, pref);
    h = testclass.make(@(x) [sin(x).*cos(x) cos(x).*exp(x)]);
    pass(n, 6) = norm(h.values - g.values, inf) < ...
        max(10*h.vscale.*h.epslevel);
    
    % Compose g(f), when f and g are CHEBTECH objects:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 7) = norm(h.values - sin(x.^2), inf) < 10*h.vscale.*h.epslevel;
    
    % Compose g(f), when f and g are CHEBTECH objects and g is array-valued:
    f = testclass.make(@(x) x.^2);
    g = testclass.make(@(x) [sin(x) cos(x)]);
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 8) = norm(h.values - [sin(x.^2) cos(x.^2)], inf) < ...
        10*max(h.vscale.*h.epslevel);
    
    % Compose g(f), when f and g are CHEBTECH objects and f is array-valued:
    f = testclass.make(@(x) [x x.^2]);
    g = testclass.make(@(x) sin(x));
    h = compose(f, g);
    x = testclass.chebpts(length(h));
    pass(n, 9) = norm(h.values - [sin(x) sin(x.^2)], inf) < ...
        10*max(h.vscale.*h.epslevel);
    
    % We cannot expect to compose two array-valued CHEBTECH objects f(g):
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) [sin(x), cos(x)]);
        compose(f, g);
        pass(n, 10) = false;
    catch ME 
        pass(n, 10) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:arrval');
    end
    
    % We cannot expect to compose two array-valued CHEBTECH objects in this way:
    try 
        f = testclass.make(@(x) [x x.^2]);
        g = testclass.make(@(x) sin(x));
        compose(f, @plus, g)
        pass(n, 11) = false;
    catch ME
        pass(n, 11) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:dim');
    end
    
    try
        f = testclass.make(@(x) sin(x));
        g = testclass.make(@(x) 100*cos(x));
        compose(g,f);
        pass(n, 12) = false;
    catch ME
        pass(n, 12) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:compose:range');
    end
end

end

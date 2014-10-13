% Test file for trigtech/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Compose a scalar-valued TRIGTECH object with sin(x):
f = testclass.make(@(x) pi*cos(pi*(x-0.1)));
g = compose(f, @sin, [], [], pref);
h = testclass.make(@(x) sin(pi*cos(pi*(x-0.1))));
pass(1) = norm(h.coeffs - g.coeffs, inf) < 10*h.vscale.*h.epslevel;

% Compose an array-valued TRIGTECH object with sin(x):
f = testclass.make(@(x) [pi*cos(pi*x) pi*cos(2*pi*x)]);
g = compose(f, @sin, [], [], pref);
h = testclass.make(@(x) [sin(pi*cos(pi*x)) sin(pi*cos(2*pi*x))]);
pass(2) = norm(h.coeffs - g.coeffs, inf) < ...
    10*max(h.vscale.*h.epslevel);

% Compose an array-valued TRIGTECH object with sin(x):
f = testclass.make(@(x) [pi*cos(pi*x) pi*cos(2*pi*x)]);
g = compose(f, @sin, [], [], pref);
x = g.points();
values = g.values;
pass(3) = norm(sin([pi*cos(pi*x) pi*cos(2*pi*x)]) - values, inf) < ...
    10*max(h.vscale.*h.epslevel);

% Compose an array-valued TRIGTECH object with sin(x):
f = testclass.make(@(x) [pi*cos(pi*x) pi*cos(2*pi*x) pi*cos(3*pi*x)]);
g = compose(f, @sin, [], [], pref);
x = g.points();
values = g.values;
pass(4) = norm(sin([pi*cos(pi*x) pi*cos(2*pi*x) pi*cos(3*pi*x)]) - values, inf) < ...
    10*max(h.vscale.*h.epslevel);

% Compose 2 TRIGTECH objects with a binary function:
f1 = testclass.make(@(x) exp(sin(pi*x)));
f2 = testclass.make(@(x) exp(cos(pi*x)));
g = compose(f1, @plus, f2, [], pref);
x = g.points;
h = testclass.make(exp(sin(pi*x)) + exp(cos(pi*x)));
hvalues = h.coeffs2vals(h.coeffs);
gvalues = g.coeffs2vals(g.coeffs);
pass(5) = norm(hvalues - gvalues, inf) < 10*h.vscale.*h.epslevel;

% Compose 2 array-valued TRIGTECH objects with a binary function:
f1 = testclass.make(@(x) exp([sin(pi*x) cos(pi*x)]));
f2 = testclass.make(@(x) exp([cos(pi*x) sin(pi*cos(pi*x))]));
g = compose(f1, @times, f2, [], pref);
x = g.points;
h = testclass.make([exp(sin(pi*x)+cos(pi*x)) exp(cos(pi*x)+sin(pi*cos(pi*x)))]);
hvalues = h.coeffs2vals(h.coeffs);
gvalues = g.coeffs2vals(g.coeffs);
pass(6) = norm(hvalues - gvalues, inf) < ...
    max(10*h.vscale.*h.epslevel);

% Compose g(f), when f and g are TRIGTECH objects:
f = testclass.make(@(x) sin(pi*x));
g = testclass.make(@(x) exp(cos(pi*x)));
h = compose(f, g); 
hvalues = h.coeffs2vals(h.coeffs);
x = testclass.trigpts(length(h));
pass(7) = norm(hvalues - exp(cos(pi*sin(pi*x))), inf) < 10*h.vscale.*h.epslevel;

% Compose g(f), when f and g are TRIGTECH objects and g is array-valued:
f = testclass.make(@(x) cos(pi*sin(pi*x)));
g = testclass.make(@(x) [sin(pi*(x-0.1)) cos(pi*(x+0.5))]);
h = compose(f, g);
x = testclass.trigpts(length(h));
hvalues = h.coeffs2vals(h.coeffs);
pass(8) = norm(hvalues - [sin(pi*(cos(pi*sin(pi*x))-0.1)) cos(pi*(cos(pi*sin(pi*x))+0.5))], inf) < ...
    10*max(h.vscale.*h.epslevel);

% Compose g(f), when f and g are TRIGTECH objects and f is array-valued:
f = testclass.make(@(x) [sin(pi*(x-0.1)) cos(pi*(x+0.5))]);
g = testclass.make(@(x) cos(pi*sin(pi*x)));
h = compose(f, g);
x = testclass.trigpts(length(h));
hvalues = h.coeffs2vals(h.coeffs);
pass(9) = norm(hvalues - cos(pi*sin(pi*[sin(pi*(x-0.1)) cos(pi*(x+0.5))])), inf) < ...
    10*max(h.vscale.*h.epslevel);

% We cannot expect to compose two array-valued TRIGTECH objects f(g):
try 
    f = testclass.make(@(x) [exp(pi*cos(pi*(x-0.14))) exp(pi*sin(pi*(x-0.14)))]);
    g = testclass.make(@(x) [sin(pi*x), cos(pi*x)]);
    compose(f, g);
    pass(10) = false;
catch ME 
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:compose:arrval');
end

% We cannot expect to compose two array-valued TRIGTECH objects in this way:
try 
    f = testclass.make(@(x) [cos(pi*x) sin(pi*(x-0.17))]);
    g = testclass.make(@(x) sin(pi*x));
    compose(f, @plus, g)
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:compose:dim');
end
    
end

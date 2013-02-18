function pass = test_flipud(pref)

if ( nargin < 1 )
    pref = funcheb2.pref;
end
tol = 10*pref.funcheb2.eps;

f = funcheb2(@(x) sin(x+.5), 0, pref);
g = funcheb2(@(x) sin(-x+.5), 0, pref);
h = flipud(f);
pass(1) = norm(g.values - h.values, inf) < tol;

f = funcheb2(@(x) [sin(x+.5), exp(x)], 0, pref);
g = funcheb2(@(x) [sin(-x+.5), exp(-x)], 0, pref);
h = flipud(f);
pass(2) = norm(g.values - h.values, inf) < tol;

f = funcheb2(@(x) sin(1i*x+.5), 0, pref);
g = funcheb2(@(x) sin(-1i*x+.5), 0, pref);
h = flipud(f);
pass(3) = norm(g.values - h.values, inf) < tol;

f = funcheb2(@(x) [sin(x+.5), exp(1i*x)], 0, pref);
g = funcheb2(@(x) [sin(-x+.5), exp(-1i*x)], 0, pref);
h = flipud(f);
pass(4) = norm(g.values - h.values, inf) < tol;

end
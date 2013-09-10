function pass = test_hypot(pref)

if ( nargin == 0 ) 
    pref = chebfun.pref();
end

% test the exmaple given in the file:

f = chebfun(@(x) 3*[1e300*x 1e-300*x], pref);
g = chebfun(@(x) 4*[1e300*x 1e-300*x], pref);
h = hypot(f, g);

seedRNG(6178);
xx = 2 * rand(100, 1) - 1;
ff = 3*[1e300*xx 1e-300*xx];
gg = 4*[1e300*xx 1e-300*xx];
hh = hypot(ff, gg);

pass(1) = norm(feval(h, xx) - hh, inf)/vscale(h) < epslevel(h);

end

function pass = test_abs(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@(x) sin(pi*x), [-3:3], pref);
g = abs(f);
h = chebfun(@(x) abs(sin(pi*x)), [-3:3], pref);
pass(1) = normest(h - g) < 10*get(f, 'epslevel')*get(f, 'hscale');

end
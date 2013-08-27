function pass = test_sign(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@(x) sin(pi*x), [-3:3], pref);
g = sign(f);
h = chebfun({-1, 1, -1, 1, -1, 1}, [-3:3]);
pass(1) = normest(h - g) < get(f, 'epslevel')*get(f, 'hscale');
p = pref;
h2 = chebfun(@(x) sign(sin(pi*x)), [-3:3], pref);
pass(2) = normest(h - g) < get(f, 'epslevel')*get(f, 'hscale');


end



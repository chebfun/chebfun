function pass = test_chebvar(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

chebvar x
f = sin(x);
g = chebfun(@sin, pref);
pass(1) = norm(f - g, inf) < epslevel(f);

chebvar x [0 10]
f = sin(x.^2) + sin(x).^2;
g = chebfun(@(x) sin(x.^2) + sin(x).^2, [0 10], pref);
pass(2) = norm(f - g, inf) < 10*epslevel(f);

chebvar x y z [0 1]
f = sin(x);
g = cos(y);
h = chebfun(@(x) sin(x) + cos(x), [0 1], pref);
pass(3) = norm((f + g) - h, inf) < epslevel(h);

end
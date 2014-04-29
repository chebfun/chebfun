% Test file for @deltafun/chebpoly.m.

function pass = test_chebpoly(pref)

if (nargin < 1)
    pref = chebpref();
end

a = -4; b = 4;

f = fun.constructor(@(x) sin(x), [a, b]);
mag = .9*(a + (b-a)*rand(3,3));
loc = .9*(a + (b-a)*rand(1,3));

df = deltafun(f, mag, loc);

pass(1) = all( chebpoly(f) == chebpoly(df) );

end
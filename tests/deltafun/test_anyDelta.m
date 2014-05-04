% Test file for @deltafun/isequal.m

function pass = test_anyDelta(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
d = deltafun();

pass(1) = ~anyDelta(d);

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = rand(1,5);
d = deltafun(f, mag, loc);
pass(2) = anyDelta(d);

d = deltafun(f,1,0);
pass(3) = anyDelta(d);

deltaTol = pref.deltaPrefs.deltaTol;
d = deltafun(f,deltaTol/2, 0);
pass(4) = ~anyDelta(d);

d = deltafun(f, [], []);
pass(5) = ~anyDelta(d);

end
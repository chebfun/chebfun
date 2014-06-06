% Test file for @deltafun/isequal.m

function pass = test_iszero(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
d1 = deltafun();
d2 = deltafun(fun.constructor(0), []);

pass(1) = iszero(d1) && iszero(d2);

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = rand(1,5);
d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
pass(2) = ~iszero(d);


end

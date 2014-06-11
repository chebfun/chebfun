% Test file for @deltafun/chebcoeffs.m.

function pass = test_chebcoeffs(pref)

if (nargin < 1)
    pref = chebfunpref();
end

a = -4; b = 4;

f = fun.constructor(@(x) sin(x), struct('domain', [a, b]));
mag = .9*(a + (b-a)*rand(3,3));
loc = .9*(a + (b-a)*rand(1,3));

df = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));

pass(1) = all( chebcoeffs(f) == chebcoeffs(df) );

end

% Test file for @deltafun/conv.m

function pass = test_conv(pref)

if (nargin < 1)
    pref = chebfunpref();
end

%%
tol = pref.deltaPrefs.deltaTol;

d = deltafun();
f = fun.constructor(@(x) x);
d1 = deltafun(f, []);
d2 = deltafun(bndfun(0), struct('deltaMag', 1, 'deltaLoc', 0));

pass(1) = isempty(conv(d, d)) && isempty(conv(d, d1)) && isempty(conv(d1, d));

%%
g = conv(d1, d2);
g = g{1};
pass(2) = max(abs(minandmax(g-f))) < tol;

%%
g = conv(d1, diff(d2));
g = g{1};
pass(3) = max(abs(minandmax(g-1))) < tol;

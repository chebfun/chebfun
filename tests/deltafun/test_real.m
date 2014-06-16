% Test file for @deltafun/real.m

function pass = test_real(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
d = deltafun();
pass(1) = isempty(real(d));

f = bndfun(@sin);
d = deltafun(f, struct('deltaMag', 1i, 'deltaLoc', 0));
pass(2) = ~isa(real(d), 'deltafun');
h = real(1i*d);
pass(3) = isa(h, 'deltafun');
pass(4) = h.deltaMag == -1;
%%

f = fun.constructor(@(x) exp(x), struct('domain', [-1, 1]));
A = rand(4,4);
B = rand(4,4);
d = deltafun(f, struct('deltaMag', A + 1i*B, 'deltaLoc', [-.5, -.25, 0, 1]));
h = real(d);
pass(5) = all(all(h.deltaMag == A));
end

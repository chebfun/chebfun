% Test file for @chebfun/firpm.m
function pass = test_firpm(pref)
%%
f = chebfun(@(x) heaviside(x) - heaviside(x-.5), [0, 1], 'splitting', 'on' );
f.pointValues = [1 1 0].';
freqs = [0, .4, .6, 1];
n = 70;
[hc, erc] = firpm(n, freqs, f);
[hm, erm] = firpm(2*n, freqs, f(freqs));
hm = chebfun(hm.', 'coeffs', 'trig' );
fc = chebfun(hc);
fc1 = restrict(fc, [freqs(1), freqs(2)]);
fc2 = restrict(fc, [freqs(3), freqs(4)]);

fm = chebfun(hm);
fm1 = restrict(fm, [freqs(1), freqs(2)]);
fm2 = restrict(fm, [freqs(3), freqs(4)]);

[norm(fc1-1, inf),  norm(fc2, inf);
 norm(fm1-1, inf),  norm(fm2, inf);]

end

% Test file for @deltafun/sum.m

function pass = test_sum(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
% Get the tolerance:
tol = pref.deltaPrefs.deltaTol;

d = deltafun();

pass(1) = sum(d) == 0; 

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = rand(1,5);

d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));

pass(2) = norm(sum(d) - (exp(1) - exp(-1) + sum(mag(1,:))), inf) < tol;

f = fun.constructor(@(x) sin(pi*x));
d = deltafun(f, struct('deltaMag', [-1, 1], 'deltaLoc', [-1, 1]));
pass(3) = norm(sum(d)-0, inf) < tol;


end

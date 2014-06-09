% Test file for @deltafun/times.m.

function pass = test_diff(pref)

if (nargin < 1)
    pref = chebfunpref();
end

tol = pref.deltaPrefs.deltaTol;
%%
a = -4;
b = 4;
d = deltafun();
pass(1) = isempty(diff(d)) && isempty(diff(d, 4));
f = bndfun(@sin, struct('domain', [a, b]));
d = deltafun(f, struct('deltaMag', 1, 'deltaLoc', 0));
pass(2) = isequal(d, diff(d,0));

f = fun.constructor(@(x) exp(x), struct('domain', [a, b]));
mag = rand(4,4);
mag(4,4) = 1; % make sure this is not a trivial row.
loc = sort(rand(1,4));
d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
dp4 = diff(d,4);
pass(3) = iszero(diff(f,4) - dp4.funPart);
pass(4) = norm(dp4.deltaLoc - loc, inf) == 0;
A = [zeros(4,4); mag] - dp4.deltaMag;
pass(5) = ~any(A(:));

end

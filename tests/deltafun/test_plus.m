% Test file for @deltafun/plus.m.

function pass = test_plus(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
a = -4; b = 4;
% Get the tolerances:
proximityTol = pref.deltaPrefs.proximityTol;
deltaTol = pref.deltaPrefs.deltaTol;

f1 = fun.constructor(@(x) sin(x), struct('domain', [a, b]));
d1 = .9*(a + (b-a)*rand(3,3));
l1 = .9*(a + (b-a)*rand(1,3));

f2 = fun.constructor(@(x) cos(x), struct('domain', [a, b]));
d2 = .9*(a + (b-a)*rand(3,3));
l2 = .9*(a + (b-a)*rand(1,3));

df1 = deltafun(f1, struct('deltaMag', d1, 'deltaLoc', l1));
df2 = deltafun(f2, struct('deltaMag', d2, 'deltaLoc', l1));

pass(1) = isempty(deltafun() + deltafun());
pass(2) = isempty(deltafun() + df1) && isempty(df1 + deltafun());

s = df1 + df2;
A = d1 + d2;
[l1, idx] = sort(l1);
A = A(:, idx);
pass(3) = norm(s.deltaMag - A, inf) == 0 && ...
    norm(s.deltaLoc - sort(l1), inf) == 0 && ...
    iszero( f1+f2 - s.funPart);

l = rand(1, 5);
d = rand(3, 5);
[sl, idx] = sort(l);
A = d(:, idx);
s = deltafun(f1, struct('deltaMag', d, 'deltaLoc', l)) ...
    + deltafun(f1, struct('deltaMag', A, 'deltaLoc', sl));
pass(4) = norm(s.deltaMag - 2*A, inf) == 0 && ...
    norm(s.deltaLoc - sl, inf) == 0 && ...
    iszero( 2*f1 - s.funPart);

data1.deltaMag = [1 2 3 4];
data1.deltaLoc = [-.25 .5 -.5 -.8];
df1 = deltafun(f1, data1);
data2.deltaMag = [1 2 3 4];
data2.deltaLoc = [-.25 .6 -.5 -.7];
df2 = deltafun(f2, data2);
s = df1 + df2;
pass(5) = norm(s.deltaLoc - [-.8 -.7 -.5 -.25 .5 .6], inf) == 0 && ...
    norm(s.deltaMag - [4 4 6 2 2 2], inf) == 0;
end

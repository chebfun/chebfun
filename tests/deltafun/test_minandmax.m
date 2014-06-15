% Test file for @deltafun/test_minandmax.m

function pass = test_minandmax(pref)

if (nargin < 1)
    pref = chebfunpref();
end

%%
tol = pref.deltaPrefs.deltaTol;

d = deltafun();
[vals, pos] = minandmax(d);
pass(1) = isempty(vals) && isempty(pos);

%%
f = fun.constructor(@(x) exp(x));
mag = [1 2 3 4 5;
       randn(1, 5)];
loc = sort(rand(1,5));

d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
[vals, pos] = minandmax(d);
pass(2) = norm(vals(1) - exp(-1), inf) < tol;
pass(3) = isinf(vals(2)) && vals(2) > 0;
pass(4) = norm(pos(1) - (-1), inf) < tol;
pass(5) = norm(pos(2) - loc(1), inf) < tol;
%%
mag = [-1 -2 -3 -4 -5;
       rand(1, 5)];
loc = sort(rand(1,5));

d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
[vals, pos] = minandmax(d);
pass(6) = isinf(vals(1)) && vals(1) < 0;
pass(7) = norm(vals(2) - exp(1), inf) < tol;
pass(8) = norm(pos(1) - loc(1), inf) < tol;
pass(9) = norm(pos(2) - 1, inf) < tol;

%%
mag = [-1 2 -3 -4 5;
       rand(1, 5)];
loc = sort(rand(1,5));

d = deltafun(f, struct('deltaMag', mag, 'deltaLoc', loc));
[vals, pos] = minandmax(d);
pass(10) = isinf(vals(1)) && vals(1) < 0;
pass(11) = isinf(vals(2)) && vals(2) > 0;
pass(12) = norm(pos(1) - loc(1), inf) < tol;
pass(13) = norm(pos(2) - loc(2), inf) < tol;

end

% Test file for @deltafun/feval.m

function pass = test_feval(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

seedRNG(1337)

%%
f = bndfun(@sin);
d = deltafun(f, struct('deltaMag', 1, 'deltaLoc', 0));
pass(1) = isinf(feval(d, 0)) && feval(d, 0) > 0;
pass(2) = isinf(feval(-d, 0)) && feval(-d, 0) <0;

%%
f = fun.constructor(@(x) sin(x));
d = deltafun(f, struct('deltaMag', [], 'deltaLoc', []));
x = rand(1, 4);
pass(2) = norm(feval(f, x) - feval(d, x), inf) == 0;

%%
x = rand(1,4);
d = deltafun(f, struct('deltaMag', rand(1, 4), 'deltaLoc', x));
pass(3) = all(isinf(feval(d, x)));

%%
x = chebfun('x');
d = dirac(x-1);
pass(4) = isinf(feval(d, 1));
val = feval(-d, 1);
pass(5) = isinf(val) && val < 0;
pass(6) = isinf(feval(d, 'right'));
pass(7) = feval(d, 'left') == 0;
pass(8) = feval(d, 1, 'left' ) == 0;
pass(9) = feval(d, 0, 'left' ) == 0;
pass(10) = feval(d, 0, 'right' ) == 0;
d = diff(heaviside(x));
pass(11) = isinf(feval(d, 0));
pass(12) = feval(d, 'right') == 0;
pass(13) = feval(d, 'left') == 0;
pass(14) = feval(d, 0, 'left' ) == 0;
pass(15) = feval(d, 0, 'right' ) == 0;

end

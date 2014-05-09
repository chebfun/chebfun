% Test file for @deltafun/feval.m

function pass = test_feval(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end
%%
f = bndfun(@sin);
d = deltafun(f, 1,0);
pass(1) = isinf(feval(d, 0)) && feval(d, 0) > 0;
pass(2) = isinf(feval(-d, 0)) && feval(-d, 0) <0;

%%
f = fun.constructor(@(x) sin(x));
d = deltafun(f, [], []);
x = rand(1, 4);
pass(2) = norm(feval(f, x) - feval(d, x), inf) == 0;

%%
x = rand(1,4);
d = deltafun(f, rand(1,4), x);
pass(3) = all(isinf(feval(d, x)));

end
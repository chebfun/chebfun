% Test file for @deltafun/cumsum.m

function pass = test_cumsum(pref)

if (nargin < 1)
    pref = chebpref();
end

%%
tol = pref.deltaPrefs.deltaTol;

d = deltafun();
F = cumsum(d);
pass(1) = isempty(F);

%%
f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = sort(rand(1,5));

d = deltafun(f, mag, loc);
[F, rval] = cumsum(d);
pass(2) = iscell(F);
%%
f = fun.constructor(@(x) sin(pi*x));
d = deltafun( f, [-1, 1], [-1, 1]);
[F, rJump] = cumsum(d);

pass(3) = ~isa(F, 'deltafun') && feval(F, -1) == -1 ...
    && rJump == 1 && feval(F, 1) == -1; 

end
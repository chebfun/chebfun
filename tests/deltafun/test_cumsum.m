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

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = sort(rand(1,5));

d = deltafun(f, mag, loc);
[F, jumpVals] = cumsum(d);

pass(2) = max(abs(F.funPart - cumsum(f))) < tol;

f = fun.constructor(@(x) sin(pi*x));
d = deltafun( f, [-1, 1], [-1, 1]);
[F, jumpVals, locations] = cumsum(d);

pass(3) = max(abs(F.funPart - cumsum(f))) < tol && ...
    norm(jumpVals - [-1, 1], inf) < tol && ...
    norm(locations - [-1, 1], inf) < tol;

%%
clc
f = fun.constructor(@(x) sin(x));
d = deltafun(f, rand(2,5), rand(1,5));
d = deltafun(1,0)
D = cumsum(d);

for i = 1:length(D)
    if ( isa( D{i}, 'deltafun') )
        plot( D{i}.funPart );
    else
        plot(D{i});
    end
    hold on
end
hold off
shg
d.deltaMag
d.location

%%
f = fun.constructor(0);
d = deltafun(f, [1 1 1 1; 1 1 1 1], [-1, -.5, .5, 1] );
D = cumsum(d);
for i = 1:length(D)
    if ( isa( D{i}, 'deltafun') )
        plot( D{i}.funPart );
        D{i}.deltaMag
        D{i}.location
    else
        plot(D{i});
    end
    hold on
end
hold off
%%
end
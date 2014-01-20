% Test file for @deltafun/cumsum.m

function pass = test_cumsum(pref)

% if (nargin < 1)
%     pref = chebpref();
% end

%%
tol = deltafun.pref.deltafun.deltaTol;

d = deltafun();
F = cumsum(d);
pass(1) = iszero(F.funPart) && isempty(F.impulses);

f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = sort(rand(1,5));

d = deltafun(f, mag, loc);
[F, jumpVals, locations] = cumsum(d);

idx = abs(mag(1, :)) > tol;

pass(2) = max(abs(F.funPart - cumsum(f))) < tol && ...
    norm(jumpVals - mag(1, idx), inf) < tol && ...
    norm(locations - loc(idx), inf) < tol;

f = fun.constructor(@(x) sin(pi*x));
d = deltafun( f, [-1, 1], [-1, 1]);
[F, jumpVals, locations] = cumsum(d);

pass(3) = max(abs(F.funPart - cumsum(f))) < tol && ...
    norm(jumpVals - [-1, 1], inf) < tol && ...
    norm(locations - [-1, 1], inf) < tol;

end
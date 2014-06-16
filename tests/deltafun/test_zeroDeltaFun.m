% Test file for @deltafun/zeroDeltaFun.m.

function pass = test_zeroDeltaFun(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
dTol = pref.deltaPrefs.deltaTol;

d = deltafun.zeroDeltaFun();
pass(1) = iszero(d.funPart);
pass(2) = all( d.funPart.domain == [-1, 1] );
pass(3) = isempty(d.deltaMag);
pass(4) = isempty(d.deltaLoc);

d = deltafun.zeroDeltaFun([4, 5]);
pass(5) = iszero(d.funPart);
pass(6) = all( d.funPart.domain == [4, 5] );
pass(7) = isempty(d.deltaMag);
pass(8) = isempty(d.deltaLoc);

end

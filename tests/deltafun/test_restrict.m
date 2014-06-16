% Test file for @deltafun/restrict.m

function pass = test_restrict(pref)

if (nargin < 1)
    pref = chebfunpref();
end
%%
d = deltafun();
pass(1) = isempty(restrict(d, [-.5, .5]));

f = bndfun(@sin);
d = deltafun(f, struct('deltaMag', 1, 'deltaLoc', 0));
pass(2) = ~isa(restrict(d, [-1, -.5]), 'deltafun');
pass(3) = ~isa(restrict(d, [.5, 1]), 'deltafun');
pass(4) = anyDelta(restrict(d, [-.5, .5]));

f = fun.constructor(@(x) exp(x), struct('domain', [-1, 1]));
d = deltafun(f, struct('deltaMag', [1 1 1 1], 'deltaLoc', [-.5, -.25, 0, 1]));
A = restrict(d, [-1 0 .5 1]);
d1 = A{1};
d2 = A{2};
d3 = A{3};

pass(5) = all(d1.deltaLoc == [-.5 -.25 0] );
pass(6) = all(d2.funPart.domain == [0, .5] );
pass(7) = all(d3.deltaMag == 1 );

%% make sure interior delta functions exactly at break points
% get divided equally in each adjacent deltafun.
f = fun.constructor(@(x) exp(x), struct('domain', [-1, 1]));
data.deltaMag = [1 1 1 1 1 1 1];
data.deltaLoc = [-1, -.5, -.25, 0, .25, .5, 1];
d = deltafun(f, data);
A = restrict(d, [-1 0 .5 1]);
d1 = A{1};
d2 = A{2};
d3 = A{3};

pass(8) = all(d1.deltaLoc == [-1 -.5 -.25 0] ) && all(d1.deltaMag == [1 1 1 .5] );
pass(9) = all(d2.deltaLoc == [0 .25 .5] ) && all(d2.deltaMag == [.5 1 .5] );
pass(10) = all(d3.deltaLoc == [.5 1] ) && all(d3.deltaMag == [.5 1] );
end

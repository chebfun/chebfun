% Test file for @deltafun/isempty.m

function pass = test_isempty(pref)

if (nargin < 1)
    pref = chebpref();
end
%%
d = deltafun();
pass(1) = isempty(d);

f = bndfun([]);
d = deltafun(f, [], []);
pass(2) = isempty(d);

d = deltafun(f, [], []);
pass(3) = isempty(d);

d = deltafun(f, [], [], []);
pass(4) = isempty(d);

end
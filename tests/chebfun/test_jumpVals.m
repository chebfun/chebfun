function pass = test_jumpVals(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@(x) x, [-1 0 .5 1]);
pass(1) = all(chebfun.jumpVals(f.funs, f.domain) == f.domain);

pass(2) = all(chebfun.jumpVals(f.funs, f.domain, {@(x) x, @(x) x, 1}) == f.domain);

pass(3) = all(chebfun.jumpVals(f.funs, f.domain, @(x) x + 100*(x==0)) == [-1 100 .5 1]);

pass(4) = all(chebfun.jumpVals(f.funs, f.domain, @sign) == [-1 0 1 1]);


end
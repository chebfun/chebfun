function pass = test_minandmax2est(pref)
% Test MINANDMAX2EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) x, [ -3, 4, -1, 2 ]);
mM = minandmax2est(f);
pass(1) = isSubset(mM, [ -3, 4], tol);

end

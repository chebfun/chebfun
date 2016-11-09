function pass = test_minandmax2est(pref)
% Test DISKFUN/MINANDMAX2EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) x);
mM = minandmax2est(f);
pass(1) = ( isSubset(mM, [-1, 1], tol) );

end
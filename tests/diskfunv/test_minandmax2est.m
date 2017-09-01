function pass = test_minandmax2est(pref)
% Test MINANDMAX2EST()

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

F = diskfunv(@(x,y) x, @(x,y) y);
rangeF = minandmax2est(F);
pass(1) = isSubset(rangeF, [-1, 1, -1, 1], tol);

end
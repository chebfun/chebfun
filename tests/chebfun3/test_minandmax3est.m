function pass = test_minandmax3est(pref)
% Test MINANDMAX3EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x, [ -2, 4, -1, 1, -1, 1 ]);
mM = minandmax3est(f);
pass(1) = ( length(mM) == 2 );
pass(2) = ( norm(mM - [ -2, 4 ]) < tol );

end

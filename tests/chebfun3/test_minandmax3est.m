function pass = test_minandmax3est(pref)
% Test MINANDMAX3EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) x, [ -2, 4, -1, 1, -1, 1 ]);
mM = minandmax3est(f);
pass(j) = ( length(mM) == 2 );
j = j + 1;
pass(j) = ( norm(mM - [ -2, 4 ]) < tol );

end

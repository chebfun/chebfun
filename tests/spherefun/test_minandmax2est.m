function pass = test_minandmax2est(pref)
% Test MINANDMAX2EST().

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1e3*pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = spherefun(@(x,y,z) z);
mM = minandmax2est(f);
pass(j) = ( norm(mM - [ -1, 1 ]) < tol );

end

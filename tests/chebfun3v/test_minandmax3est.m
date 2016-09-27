function pass = test_minandmax3est(pref)
% Test MINANDMAX3EST().

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3v();
pass(1) = isempty(minandmax3est(f));

% CHEBFUN3V (1 component)
f = chebfun3v(@(x,y,z) x, [ -2, 4, -1, 1, -1, 1 ]);
mM = minandmax3est(f);
pass(2) = ( length(mM) == 2 );
pass(3) = ( norm(mM - [ -2, 4 ]) < tol );

% CHEBFUN3V (3 components)
f = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z, [ -2, 4, 3, 17, -1, 42 ]);
mM = minandmax3est(f);
pass(4) = ( length(mM) == 6 );
pass(5) = ( norm(mM - [ -2, 4, 3, 17, -1, 42 ]) < tol );

end

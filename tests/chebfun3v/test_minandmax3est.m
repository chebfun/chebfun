function pass = test_minandmax3est(pref)
% Test MINANDMAX3EST().

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3v();
pass(j) = isempty(minandmax3est(f));
j = j + 1;

% CHEBFUN3V (1 component)
f = chebfun3v(@(x,y,z) x, [ -2, 4, -1, 1, -1, 1 ]);
mM = minandmax3est(f);
pass(j) = ( length(mM) == 2 );
j = j + 1;
pass(j) = ( norm(mM - [ -2, 4 ]) < tol );
j = j + 1;

% CHEBFUN3V (3 components)
f = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z, [ -2, 4, 3, 17, -1, 42 ]);
mM = minandmax3est(f);
pass(j) = ( length(mM) == 6 );
j = j + 1;
pass(j) = ( norm(mM - [ -2, 4, 3, 17, -1, 42 ]) < tol );

end

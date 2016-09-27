function pass = test_minandmax2est(pref)
% Test MINANDMAX2EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

% Empty CHEBFUN2V:
f = chebfun2v();
pass(1) = isempty(minandmax2est(f));

% CHEBFUN2V (1 component):
f = chebfun2v(@(x,y) x, [ -2, -1, 3, 6 ]);
box = minandmax2est(f);
pass(2) = ( length(box) == 2 );
pass(3) = ( norm(box - [ -2, -1 ]) < tol );

% CHEBFUN2V (2 components):
F = chebfun2v(@(x,y) x, @(x,y) y, [ -2, -1, 3, 6 ]);
box = minandmax2est(F);
pass(4) = ( length(box) == 4 );
pass(5) = ( norm(box - [ -2, -1, 3, 6]) < tol );

% CHEBFUN2V (3 components):
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [ -2, -1, 3, 6 ]);
box = minandmax2est(F);
pass(6) = ( length(box) == 6 );

end
function pass = test_minandmax2est(pref)
% Test MINANDMAX2EST()

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

% Empty CHEBFUN2V:
f = chebfun2v();
pass(j) = isempty(minandmax2est(f));
j = j + 1;

% CHEBFUN2V (1 component):
f = chebfun2v(@(x,y) x, [ -2, -1, 3, 6 ]);
box = minandmax2est(f);
pass(j) = ( length(box) == 2 );
j = j + 1;
pass(j) = ( norm(box - [ -2, -1 ]) < tol );
j = j + 1;

% CHEBFUN2V (2 components):
F = chebfun2v(@(x,y) x, @(x,y) y, [ -2, -1, 3, 6 ]);
box = minandmax2est(F);
pass(j) = ( length(box) == 4 );
j = j + 1;
pass(j) = ( norm(box - [ -2, -1, 3, 6]) < tol );
j = j + 1;

% CHEBFUN2V (3 components):
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y, [ -2, -1, 3, 6 ]);
box = minandmax2est(F);
pass(j) = ( length(box) == 6 );

end
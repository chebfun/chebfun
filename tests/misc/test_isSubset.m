function pass = test_isSubset(pref)

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 100 * pref.techPrefs.chebfuneps;

% Domains of chebfuns:
A = [ 0, 2 ];
B = [ 0, 2 ];
pass(1) = isSubset(A, B, tol);
B = [ 0, 1 ];
pass(2) = ~isSubset(A, B, tol);

A = [ -0.7, -0.5 ];
B = [ -1, 1];
pass(3) = isSubset(A, B, tol);

% Not exactly contained, but up to tolerance:
A = [-eps, 1];
B = [0, 1];
pass(4) = isSubset(A, B, tol);

% Domains of Chebfun2 objects:
A = [ -0.7, -0.5, 1, 2 ];
B = [ -1, 1, -1, 3 ];
pass(5) = isSubset(A, B, tol);
B = [ -1, 1, -1, 1 ];
pass(6) = ~isSubset(A, B, tol);

% Domains of Chebfun3 objects:
A = [ -1, 1, 0, 1, -1, 1 ];
B = [ -1, 1, -1, 1, -1, 1 ];
pass(7) = isSubset(A, B, tol);
A = [ -1, 1, -2, 1, -1, 1 ];
pass(8) = ~isSubset(A, B, tol);

end

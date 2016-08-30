function pass = test_isSubset(pref)

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 100 * pref.techPrefs.chebfuneps;
j = 1;

% Domains of chebfuns:
A = [ 0, 2 ];
B = [ 0, 2 ];
pass(j) = isSubset(A, B, tol);
j = j + 1;
B = [ 0, 1 ];
pass(j) = ~isSubset(A, B, tol);
j = j + 1;

A = [ -0.7, -0.5 ];
B = [ -1, 1];
pass(j) = isSubset(A, B, tol);
j = j + 1;

% Not exactly contained, but up to tolerance:
A = [-eps, 1];
B = [0, 1];
pass(j) = isSubset(A, B, tol);
j = j + 1;

% Domains of Chebfun2 objects:
A = [ -0.7, -0.5, 1, 2 ];
B = [ -1, 1, -1, 3 ];
pass(j) = isSubset(A, B, tol);
j = j + 1;
B = [ -1, 1, -1, 1 ];
pass(j) = ~isSubset(A, B, tol);
j = j + 1;

% Domains of Chebfun3 objects:
A = [ -1, 1, 0, 1, -1, 1 ];
B = [ -1, 1, -1, 1, -1, 1 ];
pass(j) = isSubset(A, B, tol);
j = j + 1;
A = [ -1, 1, -2, 1, -1, 1 ];
pass(j) = ~isSubset(A, B, tol);

end

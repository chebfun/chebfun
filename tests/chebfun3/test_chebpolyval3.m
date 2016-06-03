function pass = test_chebpolyval3()
% Check the chebpolyval3 command in CHEBFUN3

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.cheb3Prefs.chebfun3eps;

% check chebpolyval3 command:
seedRNG(0);
exactVals = rand(4,4,4);
constructorVals = chebpolyval3(chebfun3(exactVals));
pass(1) = norm(exactVals(:) - constructorVals(:)) < 10*tol;

% Check degrees are as expected:
A = rand(3, 4, 5);
f = chebfun3(A);
B = chebpolyval3(f);
pass(2) = norm(A(:) - B(:)) < tol;
end
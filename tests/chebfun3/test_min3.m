function pass = test_min3(pref)
% Test @chebfun3/min3 command.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

% Check the MIN3 function for Wagon's function from cheb.gallery3.
% This is a variation of Problem 4 in the SIAM 100-Digit Challenge
f = cheb.gallery3('wagon');

% Solution from the book: 
exactVal = -3.3283383456632;
exactLoc = [-0.158036820468905, 0.291023048609152, -0.289297798732570];
[val, loc] = min3(f);
pass(1) = norm(loc - exactLoc) < tol;
pass(2) = abs(val - exactVal) < tol;

end
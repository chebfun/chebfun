function pass = test_sum3(pref)
% Test file for @chebfun3t/sum3.

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e1 * pref.cheb3Prefs.chebfun3eps;

% Check empty case
pass(1) = isempty(sum3(chebfun3t()));

% Check constant
f = chebfun3t(@(x,y,z) 1);
pass(2) = abs(sum3(f) - 8) < tol;

% Runge function
f = chebfun3t(@(x,y,z) 1./(1+x.^2+y.^2+z.^2));
pass(3) = abs(sum3(f) - 4.28685406230184188268) < tol;

% Different domains
f = chebfun3t(@(x,y,z) x, [0, 1, 0, 2, 0, 3]);
pass(4) = abs(sum3(f) - 3) < tol;

f = chebfun3t(@(x,y,z) y, [0, 1, 0, 2, 0, 3]);
pass(5) = abs(sum3(f) - 6) < tol;

f = chebfun3t(@(x,y,z) z, [0, 1, 0, 2, 0, 3]);
pass(6) = abs(sum3(f) - 9) < tol;

end
function pass = test_min2(pref)
% Test the chebfun3/min2 command.

if ( nargin < 1 )
    pref = chebfunpref; 
end
tol = 1e12*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2 + y.^2 + z.^2);
g = chebfun(@(z) z.^2);

h1 = min2(f);
h2 = min2(f, []);
h3 = min2(f, [], [1 2]);
h4 = min2(f, [], [2 1]);

h5 = min2(f, [], [1 3]);
h6 = min2(f, [], [2 3]);

pass(1) = norm(h1 - g) < tol;
pass(2) = norm(h2 - g) < tol;
pass(3) = norm(h3 - g) < tol;
pass(4) = norm(h4 - g) < tol;
pass(5) = norm(h5 - g) < tol;
pass(6) = norm(h6 - g) < tol;

end
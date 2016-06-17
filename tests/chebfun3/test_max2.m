function pass = test_max2(pref)
% Test the chebfun3/max2 command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e7*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z));
g = chebfun(@(x) 1 + 0*x); 

h1 = max2(f);
h2 = max2(f, []);
h3 = max2(f, [], [1 2]);
h4 = max2(f, [], [2 1]);

h5 = max2(f, [], [1 3]);
h6 = max2(f, [], [3 1]);

h7 = max2(f, [], [2 3]);
h8 = max2(f, [], [3 2]);

pass(1) = norm(h1 - g) < tol;
pass(2) = norm(h2 - g) < tol;
pass(3) = norm(h3 - g) < tol;
pass(4) = norm(h4 - g) < tol;
pass(5) = norm(h5 - g) < tol;
pass(6) = norm(h6 - g) < tol;
pass(7) = norm(h7 - g) < tol;
pass(8) = norm(h8 - g) < tol;

end
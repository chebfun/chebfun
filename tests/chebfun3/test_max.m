function pass = test_max(pref)
% Test the chebfun3/max command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e11*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z));
g = chebfun2(@(x,y) 1 + 0*x); 

h1 = max(f); 
h2 = max(f, []); 
h3 = max(f, [], 2);
h4 = max(f, [], 3);  
h5 = max(f, [], 4);

pass(1) = norm(h1 - g) < tol;
pass(2) = norm(h2 - g) < tol;
pass(3) = norm(h3 - g) < tol;
pass(4) = norm(h4 - g) < tol;
pass(5) = norm(h5 - f) < tol;

end
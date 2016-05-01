function pass = test_min( pref ) 
% Test the chebfun3/min command. 

if ( nargin < 1 ) 
    pref = chebfunpref;
end
tol = 1e11*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2+y.^2+z.^2);
g = chebfun2(@(y,z) y.^2+z.^2);

h1 = min(f); 
h2 = min(f, []); 
h3 = min(f, [], 2);
h4 = min(f, [], 3);  
h5 = min(f, [], 4);

pass(1) = norm(h1 - g) < tol;
pass(2) = norm(h2 - g) < tol;
pass(3) = norm(h3 - g) < tol;
pass(4) = norm(h4 - g) < tol;
pass(5) = norm(h5 - f) < tol;

end
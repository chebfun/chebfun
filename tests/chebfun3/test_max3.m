function pass = test_max3(pref)
% Test @chebfun3/max3 command.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Check the MAX3 function for a cosine function
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun(@(x,y,z) 1 + 0*x); 
h1 = max3(f); 
pass(1) = norm(h1-g) < tol;

% Check the MAX function for a sine function
f = chebfun3(@(x,y,z) sin(x+y+z)); 
h2 = max3(f);
pass(2) = norm(h2-g) < tol;

end
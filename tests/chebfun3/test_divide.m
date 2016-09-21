function pass = test_divide(pref)
% Test contour

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun3(@(x,y,z) cos(x.*y.*z)./2); 

pass(1) = norm( f./2 - g ) < tol;

pass(2) = norm( f/2 - g ) < tol;

pass(3) = norm( 2.\f - g ) < tol;

pass(4) = norm( 2\f - g ) < tol;

end
function pass = test_divide(pref)
% Test contour

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun3(@(x,y,z) cos(x.*y.*z)./2); 

pass(j) = norm( f./2 - g ) < tol; j = j + 1; 
pass(j) = norm( f/2 - g ) < tol; j = j + 1; 
pass(j) = norm( 2.\f - g ) < tol; j = j + 1; 
pass(j) = norm( 2\f - g ) < tol; j = j + 1; 

end
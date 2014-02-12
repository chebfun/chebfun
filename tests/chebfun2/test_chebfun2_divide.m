function pass = test_chebfun2_divide( pref ) 
% Test contour

if ( nargin == 0) 
    pref = chebpref; 
end

tol = 1000*pref.cheb2Prefs.eps; 
j = 1; 

f = chebfun2(@(x,y) cos(x.*y)); 
g = chebfun2(@(x,y) cos(x.*y)./2); 

pass(j) = norm( f./2 - g ) < tol; j = j + 1; 
pass(j) = norm( f/2 - g ) < tol; j = j + 1; 
pass(j) = norm( 2.\f - g ) < tol; j = j + 1; 
pass(j) = norm( 2\f - g ) < tol; j = j + 1; 


end
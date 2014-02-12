function pass = test_chebfun2_transpose( pref ) 
% Test transpose

if ( nargin == 0) 
    pref = chebpref; 
end

tol = 1000*pref.cheb2Prefs.eps; 
j = 1; 

% symmetric function:
f = chebfun2(@(x,y) cos(x.*y)); 
pass(j) = norm( f - f') < tol; j = j + 1; 
pass(j) = norm( f - f.') < tol; j = j + 1; 

f = chebfun2(@(x,y) cos(x.*y) + x, [-3 4 -1 0]); 
g = chebfun2(@(x,y) cos(x.*y) + y, [-1 0 -3 4]); 
pass(j) = norm( f' - g ) < tol; j = j + 1;
pass(j) = norm( f.' - g ) < tol; j = j + 1;

f = chebfun2(@(x,y) cos(x.*y) + 1i*x, [-3 4 -1 0]); 
g1 = chebfun2(@(x,y) cos(x.*y) - 1i*y, [-1 0 -3 4]); 
g2 = chebfun2(@(x,y) cos(x.*y) + 1i*y, [-1 0 -3 4]); 
pass(j) = norm( f' - g1 ) < tol; j = j + 1;
pass(j) = norm( f.' - g2) < tol; j = j + 1;

end
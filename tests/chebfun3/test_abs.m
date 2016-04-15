function pass = test_abs(pref)
% Test abs. 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z) + 2); 
pass(j) = norm(f - abs(f)) < tol; j = j + 1; 
pass(j) = norm(f - abs(-f)) < tol; j = j + 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z) + 2, [-3 4 -1 1 -2 0]); 
pass(j) = norm(f - abs(f)) < tol; j = j + 1; 
pass(j) = norm(f - abs(-f)) < tol; j = j + 1; 

end
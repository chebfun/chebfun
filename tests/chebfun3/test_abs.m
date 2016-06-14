function pass = test_abs(pref)
%Test chebfun3/abs.

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps; 

f = chebfun3(@(x,y,z) cos(x.*y.*z) + 2); 
pass(1) = norm(f - abs(f)) < tol;
pass(2) = norm(f - abs(-f)) < tol;

f = chebfun3(@(x,y,z) cos(x.*y.*z) + 2, [-3 4 -1 1 -2 0]); 
pass(3) = norm(f - abs(f)) < tol;
pass(4) = norm(f - abs(-f)) < tol;

try
    f = chebfun3(@(x,y,z) .9*x + y.^2);
    g = abs(f);
catch
    pass(5) = 1;
end

end
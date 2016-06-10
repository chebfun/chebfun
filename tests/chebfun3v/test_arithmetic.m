function pass = test_arithmetic(pref)
% Check the Chebfun3v constructor for simple arithmetic operations.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;
j = 1;

f = @(x,y,z) cos(x.*y.*z); 
f=chebfun3v(f,f);

g = @(x,y,z) sin(y); 
g=chebfun3v(g,g);

% exact answers:
plus_exact = @(x,y,z) cos(x.*y.*z) + sin(y);
plus_exact = chebfun3v(plus_exact, plus_exact);

minus_exact = @(x,y,z) cos(x.*y.*z) - sin(y); 
minus_exact = chebfun3v(minus_exact, minus_exact);

mult_exact = @(x,y,z) cos(x.*y.*z).*sin(y); 
mult_exact=chebfun3v(mult_exact, mult_exact);

pass(j) = norm(f + g - plus_exact) < tol; 
j=j+1;

pass(j) = norm(f - g - minus_exact) < tol; 
j=j+1;

pass(j) = norm(f.*g - mult_exact) < tol; 
j=j+1;

end
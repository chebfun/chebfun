function pass = test_chebfun2v_arithmetic
% Check the Chebfun2v constructor for simple arithmetic operations.
% Alex Townsend, March 2013.

% These function chosen so that scl does not change.
f = @(x,y) cos(x); f=chebfun2v(f,f);
g = @(x,y) sin(y); g=chebfun2v(g,g);
% exact answers.
plus_exact = @(x,y) cos(x) + sin(y); 
plus_exact=chebfun2v(plus_exact, plus_exact);
minus_exact = @(x,y) cos(x) - sin(y); 
minus_exact=chebfun2v(minus_exact, minus_exact);
mult_exact = @(x,y) cos(x).*sin(y); 
mult_exact=chebfun2v(mult_exact, mult_exact);
pow_exact = @(x,y) cos(x).^sin(y); 
pow_exact=chebfun2v(pow_exact, pow_exact);

tol = 1e-14;

pass(1) = norm(f + g - plus_exact) < tol;
pass(2) = norm(f - g - minus_exact) < tol;
pass(3) = norm(f.*g - mult_exact) < tol;
pass(4) = norm(f.^g - pow_exact) < tol;
pass = all(pass);

end

function pass = test_arithmetic( pref )
% Check the Chebfun2v constructor for simple arithmetic operations.
% Alex Townsend, March 2013.


if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

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

pass(j) = norm(f + g - plus_exact) < tol; j=j+1;
pass(j) = norm(f - g - minus_exact) < tol; j=j+1;
pass(j) = norm(f.*g - mult_exact) < tol; j=j+1;

end

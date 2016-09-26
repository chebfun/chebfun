function pass = test_compose(pref)
% Check that CHEBFUN3T composition operations are working.

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3t(@(x,y,z) cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 

% Multiplication
exact = @(x,y,z) (cos(x.*y.*z) + sin(x.*y.*z) + y -.1) .* ...
    sin((x-.1).*(y+.4).*(z+.8));
g = chebfun3t(@(x,y,z) exact(x,y,z)); 
pass(1) = norm(g - chebfun3t(@(x,y,z) f(x,y,z).*sin((x-.1).*(y+.4).*(z+.8)))) < tol; 

% Sine
exact = @(x,y,z) sin(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3t(@(x,y,z) exact(x,y,z));
pass(2) = norm(g - sin(f)) < tol; 

% Cosine
exact = @(x,y,z) cos(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3t(@(x,y,z) exact(x,y,z)); 
pass(3) = norm(g - cos(f)) < tol;

% Sinh
exact = @(x,y,z) sinh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3t(@(x,y,z) exact(x,y,z));
pass(4) = norm(g - sinh(f)) < tol;

% Cosh
exact = @(x,y,z) cosh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3t(@(x,y,z) exact(x,y,z)); 
pass(5) = norm( g - cosh(f) ) < tol;

% Tanh
exact = @(x,y,z) tanh(-(cos(x.*y.*z) + sin(x.*y.*z) + y -.1)); 
g = chebfun3t(@(x,y,z) exact(x,y,z)); 
pass(6) = norm(g - tanh(-f)) < tol;

% Exp
exact = @(x,y,z) exp(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3t(@(x,y,z) exact(x,y,z)); 
pass(7) = norm(g - exp(f)) < tol;

% Multiple operations: 
f = chebfun3t(@(x,y,z) sin(10*x.*y.*z), [-1 2 -1 1 -3 -1]);
pass(8) = norm(f+f+f-3*f) < 100*tol; 

pass(9) = norm(f.*f-f.^2) < tol;

end
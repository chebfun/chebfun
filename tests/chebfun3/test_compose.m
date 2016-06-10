function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 

% Multiplication
exact = @(x,y,z) (cos(x.*y.*z) + sin(x.*y.*z) + y -.1) .* ...
    sin((x-.1).*(y+.4).*(z+.8));
g = chebfun3(@(x,y,z) exact(x,y,z)); 
[x, y, z] = cheb.xyz;
pass(j) = norm(g - f.*sin((x-.1).*(y+.4).*(z+.8))) < tol; 
j = j + 1;

% Sine
exact = @(x,y,z) sin(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z));
pass(j) = norm(g - sin(f)) < tol; 
j = j + 1;

% Cosine
exact = @(x,y,z) cos(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(j) = norm(g - cos(f)) < tol; 
j = j + 1;

% Sinh
exact = @(x,y,z) sinh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3(@(x,y,z) exact(x,y,z));
pass(j) = norm(g - sinh(f)) < tol; 
j = j + 1;

% Cosh
exact = @(x,y,z) cosh(cos(x.*y.*z) + sin(x.*y.*z) + y -.1);
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(j) = norm( g - cosh(f) ) < tol; 
j = j + 1;

% Tanh
exact = @(x,y,z) tanh(-(cos(x.*y.*z) + sin(x.*y.*z) + y -.1)); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(j) = norm(g - tanh(-f)) < tol; 
j = j + 1;

% Exp
exact = @(x,y,z) exp(cos(x.*y.*z) + sin(x.*y.*z) + y -.1); 
g = chebfun3(@(x,y,z) exact(x,y,z)); 
pass(j) = norm(g - exp(f)) < tol; 
j = j + 1;

% Multiple operations: 
f = chebfun3(@(x,y,z) sin(10*x.*y.*z), [-1 2 -1 1 -3 -1]);
pass(j) = norm(f+f+f-3*f) < 100*tol;
j=j+1; 

pass(j) = norm(f.*f-f.^2) < tol;

end
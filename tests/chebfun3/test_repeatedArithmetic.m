function pass = test_repeatedArithmetic(pref)

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

% Add 50 functions together: 
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun3(0);
for jj = 1:50 
    g = g + f; 
end
pass(1) = norm(g - 50*f ) < 10 * tol; 

% Multiply 10 functions together: 
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = f;
for jj = 1:10 
    g = g .* f; 
end
pass(2) = norm(g - f.^11) < tol;

% Multiply 20 functions together: 
f = chebfun3(@(x,y,z) sin(x.*y.*z)); 
g = f;
for jj = 1:20 
    g = g .* f; 
end
pass(3) = norm(g - f.^21) < tol;

end
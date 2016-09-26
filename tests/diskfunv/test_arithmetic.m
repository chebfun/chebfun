function pass = test_arithmetic( pref )
% Check simple arithmetic operations in diskfunv

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3*pref.cheb2Prefs.chebfun2eps;
j = 1;

f = @(x,y) cos(2*pi*x.*y);
g = @(x,y) sin(2*pi*x.*y);
h = @(x,y) 1 + x.*y;
u = diskfunv(f,g);
v = diskfunv(h,f);

% Exact answers.
plus_fh = @(x,y) f(x,y) + h(x,y); 
plus_gf = @(x,y) g(x,y) + f(x,y);  
plus_exact = diskfunv( plus_fh, plus_gf);

minus_fh = @(x,y) f(x,y) - h(x,y); 
minus_gf = @(x,y) g(x,y) - f(x,y); 
minus_exact = diskfunv(minus_fh, minus_gf);

mult_fh = @(x,y) f(x,y).*h(x,y); 
mult_gf = @(x,y) g(x,y).*f(x,y); 
mult_exact = diskfunv(mult_fh, mult_gf);

pass(j) = norm(u + v - plus_exact) < tol; j=j+1;
pass(j) = norm(u - v - minus_exact) < tol; j=j+1;
pass(j) = norm(u.*v - mult_exact) < tol; 

end

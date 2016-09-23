function pass = test_arithmetic( pref )
% Check simple arithmetic operations in spherefunv

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3*pref.cheb2Prefs.chebfun2eps;
j = 1;

f = @(x,y,z) cos(2*pi*x.*y.*z);
g = @(x,y,z) sin(2*pi*x.*y.*z);
h = @(x,y,z) 1 + x.*y.*z;
u = spherefunv(f,g,h);
v = spherefunv(h,f,g);

% Exact answers.
plus_fh = @(x,y,z) f(x,y,z) + h(x,y,z); 
plus_gf = @(x,y,z) g(x,y,z) + f(x,y,z); 
plus_hg = @(x,y,z) h(x,y,z) + g(x,y,z); 
plus_exact = spherefunv(plus_fh, plus_gf, plus_hg);

minus_fh = @(x,y,z) f(x,y,z) - h(x,y,z); 
minus_gf = @(x,y,z) g(x,y,z) - f(x,y,z); 
minus_hg = @(x,y,z) h(x,y,z) - g(x,y,z); 
minus_exact = spherefunv(minus_fh, minus_gf, minus_hg);

mult_fh = @(x,y,z) f(x,y,z).*h(x,y,z); 
mult_gf = @(x,y,z) g(x,y,z).*f(x,y,z); 
mult_hg = @(x,y,z) h(x,y,z).*g(x,y,z); 
mult_exact = spherefunv(mult_fh, mult_gf, mult_hg);

pass(j) = norm(u + v - plus_exact) < tol; j=j+1;
pass(j) = norm(u - v - minus_exact) < tol; j=j+1;
pass(j) = norm(u.*v - mult_exact) < tol;

end

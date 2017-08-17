function pass = test_poisson( pref ) 
% Test Poisson solver in Chebfun2 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1e6*pref.cheb2Prefs.chebfun2eps;

% Compare against chebop2: 
m = 101; n = 101; 
N = chebop2( @(u) lap(u) ); 
N.bc = 0; 
f = chebfun2(@(x,y) 1+0*x); 
u = chebfun2.poisson(f, 0, m, n);
v = solvepde(N, f, m, n);
pass(1) = ( norm(u-v) < tol );

% Rectangular discretization sizes:
m = 87; n = 101; 
N = chebop2( @(u) lap(u) ); 
N.bc = 0; 
f = chebfun2(@(x,y) 1+0*x); 
u = chebfun2.poisson(f, 0, m, n);
v = solvepde(N, f, m, n);
pass(2) = ( norm(u-v) < tol );

% Nonstandard rectangular domain:  
m = 99; n = 101; 
a = -3; b = 1; c = 4; d = 4.2; 
N = chebop2( @(u) lap(u), [a b c d] ); 
N.bc = 0; 
f = chebfun2(@(x,y) 1+0*x, [a b c d]); 
u = chebfun2.poisson(f, 0, m, n);
v = solvepde(N, f, m, n);
pass(3) = ( norm(u-v) < tol );

% Nonzero, but constant, Dirichlet data:
m = 99; n = 101; 
a = -3; b = 1; c = 4; d = 4.2; 
N = chebop2( @(u) lap(u), [a b c d] ); 
N.bc = 1; 
f = chebfun2(@(x,y) 1+0*x, [a b c d]); 
u = chebfun2.poisson(f, 1, m, n);
v = solvepde(N, f, m, n);
pass(4) = ( norm(u-v) < tol );

% General Dirichlet data, given as function handle:
m = 201; n = 201; 
a = -3; b = 1; c = 4; d = 4.2; 
p = @(x,y) x.*y + cos(3*x.^2.*(y-.2));
N = chebop2( @(u) lap(u), [a b c d] ); 
N.lbc = @(y) p(a,y);
N.rbc = @(y) p(b,y);
N.dbc = @(x) p(x,c);
N.ubc = @(x) p(x,d);
f = chebfun2(@(x,y) 1+0*x, [a b c d]); 
u = chebfun2.poisson(f, p, m, n);
v = solvepde(N, f, m, n);
pass(5) = ( norm(u-v) < 100*tol );

% General Dirichlet data, given as chebfun2:
m = 201; n = 201; 
a = -3; b = 1; c = 4; d = 4.2; 
p = @(x,y) x.*y + cos(3*x.^2.*(y-.2));
N = chebop2( @(u) lap(u), [a b c d] ); 
N.lbc = @(y) p(a,y);
N.rbc = @(y) p(b,y);
N.dbc = @(x) p(x,c);
N.ubc = @(x) p(x,d);
g = chebfun2( p, [a b c d]);
f = chebfun2(@(x,y) 1+0*x, [a b c d]); 
u = chebfun2.poisson(f, g, m, n);
v = solvepde(N, f, m, n);
pass(6) = ( norm(u-v) < 100*tol );

if ( all(pass) ) 
    pass = 1; 
end

end



function pass = test_Poisson( ) 
% Check correctness of Poisson solver on the sphere: 

% Discretization sizes: 
m = 40; 
n = 40;
tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;
 
% Example 1: 
f = spherefun(@(lam,th) -6*cos(lam).*cos(th).*sin(th)); 
exact = spherefun(@(lam,th) sin(th).*cos(th).*cos(lam));
u = spherefun.poisson(f, 0, m, n);
pass(1) = ( norm(u - exact, inf) < tol ); 

% Example 2: 
f = spherefun(@(lam,th) -4*(3*cos(th)+5*cos(3*th)).*sin(lam).*sin(th)); 
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2);
u = spherefun.poisson(f, 0, m, n);
pass(2) = ( norm(u - exact, inf) < tol ); 

% Inline example: 
f = spherefun(@(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th));
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
            sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th));
u = spherefun.poisson(f, 0, m, n);
pass(3) = ( norm(u - exact, inf) < tol );

% Check that the mean is zero.
pass(4) = ( abs(mean2(u)) < tol );

% Check that the code properly deals with a right hand side without a 
% mean of zero.
f = spherefun(@(x,y,z) 1 + x);
warning('off','CHEBFUN:SPHEREFUN:POISSON:meanRHS');
u = spherefun.poisson(f, 0, 10);
warning('on','CHEBFUN:SPHEREFUN:POISSON:meanRHS');
% This f - mean2(f)
g = spherefun(@(x,y,z) x);
% Solution with the mean of f set to zero
v = spherefun.poisson(g, 0, 10);
pass(5) = ( norm(u - v) < tol );

% Check that the code allows the mean of the solution to be set to
% something other than zero.
f = spherefun(@(x,y,z) x.*y.*z );
u = spherefun.poisson(f, 1, 10);
pass(6) = ( abs(mean2(u)-1) < tol );


end
function pass = test_Poisson( ) 
% Check correctness of Poisson solver on the sphere: 

% Discretization sizes: 
m = 40; 
n = 40;
tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;
 
% Example 1: 
f = spherefun(@(lam,th) -6*cos(lam).*cos(th).*sin(th)); 
exact = spherefun(@(lam,th) sin(th).*cos(th).*cos(lam));
u = spherefun.poisson(f, m, n);
pass(1) = ( norm(u - exact, inf) < tol ); 

% Example 2: 
f = spherefun(@(lam,th) -4*(3*cos(th)+5*cos(3*th)).*sin(lam).*sin(th)); 
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2);
u = spherefun.poisson(f, m, n);
pass(2) = ( norm(u - exact, inf) < tol ); 

% Inline example: 
f = spherefun(@(lam,th) -6*(-1+5*cos(2*th)).*sin(lam).*sin(2*th));
exact = spherefun(@(lam,th) -2*sin(lam).*sin(2*th).*sin(th).^2 -...
            sin(lam).*sin(th).*cos(th) + .5*sin(lam).*sin(2*th).*cos(2*th));
u = spherefun.poisson(f, m, n);
pass(3) = ( norm(u - exact, inf) < tol );

end
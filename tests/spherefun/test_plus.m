function pass = test_plus( ) 
% Test spherefun plus() command 

tol = 1e3*chebfunpref().techPrefs.eps;

f1 = @(x,y,z) sin(pi*x.*y);  % Strictly even/pi-periodic
f2 = @(x,y,z) sin(pi*x.*z);  % Strictly odd/anti-periodic
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = redefine_function_handle( @(x,y,z) f1(x,y,z) + f2(x,y,z) );
lambda = rand; theta = rand; 
pass(1) = abs( feval(gplus, theta, lambda) - fplus(theta, lambda) ) < tol; 

f1 = @(lam,th) exp(cos(lam-1).*sin(th).*cos(th));  % Mixed symmetric terms
f2 = @(lam,th) exp(sin(lam-0.35).*sin(th).*cos(th));  % Mixed symmetric terms
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = @(lam,th) f1(lam,th) + f2(lam,th) ;
lambda = rand; theta = rand; 
pass(2) = abs( feval(gplus, theta, lambda) - fplus(theta, lambda) ) < tol; 

% Check that compression is working: 
f = spherefun(@(x,y,z) x.^2 + y.^2 + z.^2); 
r = rank( f ); 
g = f; 
for k = 1:10
    g = g + f; 
end 
pass(3) = norm( g - 11*f ) < vscale(g)*tol; 
pass(4) = ( rank( g ) - r ) == 0; 

% Check what happens with cancellation errors: 
f = spherefun(@(x,y,z) sin(x.*y.*z)); 
g = 2*f; 
pass(5) = ( norm( g - f - f ) < tol ); 

end
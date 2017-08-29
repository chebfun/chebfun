function pass = test_jacobian( ) 
% Test jacobian

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check jacobian of an empty diskfunv is an empty diskfun.
u = diskfunv;
f = jacobian(u);
pass(1) = isempty(f) & isa(f,'diskfun');

% Check definition: 
F = diskfunv(@(x,y) cos(pi*x.*y), @(x,y) sin(pi*y).*cos(pi*x)); 
Fx = diffx(F); 
Fy = diffy(F); 
pass(2) = ( norm(jacobian(F) - (Fx(1).*Fy(2)-Fy(1).*Fx(2)) ) < tol );

%check against true value
true = diskfun(@(x,y)-pi^2*y.*sin(pi*x.*y).*cos(pi*y).*cos(pi.*x)-...
    pi^2*x.*sin(pi*x.*y).*sin(pi*y).*sin(pi*x)); 
pass(3) = ( norm(jacobian(F)-true) < 4*tol ) ; 


end
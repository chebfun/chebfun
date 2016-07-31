function pass = test_curl( ) 
% Test curl

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Check curl of an empty diskfunv is an empty diskfunv.
u = diskfunv;
f = curl(u);
pass(1) = isempty(f) & isa(f,'diskfunv');

% Check curl of a diskfunv is a diskfun.
u = diskfunv(@(t,r) r.*cos(t), @(t,r) r.*sin(t), 'polar'); 
f = curl(u);
pass(2) = isa(f,'diskfun');

% Check definition: 
F = diskfunv(@(x,y) cos(x), @(x,y) sin(y));  
curlF = diff(F(2), 1,1) - diff(F(1), 2, 1);
pass(3) = ( norm(curlF - curl(F)) < tol );

% Check curl of the zero field is zero
f = diskfun(@(x,y) 0*x);
u = diskfunv(f,f);
pass(4) = norm(curl(u)) < tol; 


% Check curl of gradient field should be zero
f = diskfun(@(x,y) cos(x+y.^2) + sin(y) + y); 
F = gradient( f ); 
pass(5) = ( norm( curl( F ) ) < 1e3*tol ); 
end
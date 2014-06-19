function pass = test_gradys_function1( pref )
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.eps; 

% Grady's function that failed: 
g = @(x,y) exp(-1./max(1 - ((x-0.02).^2 + (y-0.033).^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,101));
err = g(xx,yy)-f(xx,yy);
pass(1) = ( norm(err(:),inf ) < 2e3*tol );

end
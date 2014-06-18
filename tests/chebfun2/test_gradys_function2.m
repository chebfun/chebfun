function pass = test_gradys_function2( pref )
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.eps; 

% Another variant on Grady's function: 
g = @(x,y) exp(-((x-0.2).^2+(y-0.33).^2)./max(1 - ((x-0.2).^2 + (y-0.33).^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,101));
err = g(xx,yy)-f(xx,yy);
pass(1) = ( norm(err(:),inf ) < 2e3*tol );

end
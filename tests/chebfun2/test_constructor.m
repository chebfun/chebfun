function pass = test_constructor( pref ) 
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.eps; 

% Can we make a chebfun2: 
f = @(x,y) cos( x ) + sin( x .* y );  % simple function. 
fstr = 'cos(x) + sin(x.*y)'; % string version.

% construct with a domain.
f1 = chebfun2( f ); 
f2 = chebfun2( f, [-1 1 -1 1] );
pass(1) = ( norm(f1 - f2) < tol ); 

% construct from a string. 
f3 = chebfun2( fstr ); 
f4 = chebfun2( fstr, [-1 1 -1 1]); 
pass(2) = ( norm( f1 - f3 ) < tol ); 
pass(3) = ( norm( f3 - f4 ) < tol );

g = @(x,y) cos(x).*y + x.*sin(y); 
f = chebfun2( g );
pass(4) = ( length(f) == 2 ); 
pass(5) = ( abs( feval(f, .1, pi/6) - feval(g, .1, pi/6) ) < tol );

% George's function that failed:
f = @(x,y) 1./(1+25*x.^2.*y.^2);
ffch = chebfun2(@(x,y) f(x,y), [-2 2 -2 2]);
xx = linspace(-2,2); 
[XX,YY] = meshgrid(xx,xx);
pass(6) = ( max(max( abs(f(XX,YY) - ffch(XX,YY) ))) < 2e4*tol );

% Grady's function that failed: 

g = @(x,y) exp(-1./max(1 - ((x-0.02).^2 + (y-0.033).^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,101));
err = g(xx,yy)-f(xx,yy);
pass(7) = ( norm(err(:),inf ) < 2e3*tol );

% Another variant on Grady's function: 
g = @(x,y) exp(-((x-0.2).^2+(y-0.33).^2)./max(1 - ((x-0.2).^2 + (y-0.33).^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,101));
err = g(xx,yy)-f(xx,yy);
pass(8) = ( norm(err(:),inf ) < 2e3*tol );

% Test building Chebfun2 objects from sample data: 
seedRNG(0);
r = rand(3);
pass(9) = norm( r - chebpoly2(chebfun2(r, 'coeffs')) ) < 10*tol; 
r = rand(4);
pass(10) = norm( r - chebpoly2(chebfun2(r, 'coeffs')) ) < 10*tol;
r = rand(4);
pass(11) = norm( r - chebpolyval2(chebfun2(r)) ) < 10*tol;

g = @(x,y) exp(-((x+pi).^2+y.^2)./max(1 - ((x+pi).^2 + y.^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,1001));
err = g(xx,yy)-f(xx,yy);
pass(12) = ( norm(err(:),inf ) < 2e7*tol );

end

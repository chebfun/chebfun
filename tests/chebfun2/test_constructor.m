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

% Test that #928 is fixed.
f = chebfun2(@(x,y) 1, 'vectorize');
g = chebfun2(1);
pass(7) = norm(f - g) < tol;
f = chebfun2(@(x,y) 1, 'vectorize');
g = chebfun2(1);
pass(8) = norm(f - g) < tol;

% Test that the 2nd argument can be a preference:
p = pref;
p.tech = @trigtech;
f = chebfun2(@(x,y) sin(pi*x).*cos(pi*y), p);
g = chebfun2(@(x,y) sin(pi*x).*cos(pi*y), [-1 1 -1 1], p);
pass(9) = norm(f - g) < tol;

end

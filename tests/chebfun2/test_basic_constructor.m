function pass = test_basic_constructor( pref ) 
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb2Prefs.eps; 

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
pass(6) = ( max(max( abs(f(XX,YY) - ffch(XX,YY) ))) < 2e3*tol );

%%
% Constructing from matrix of values 
% A = randn( 100 ); 
% f = chebfun2( A ); 
% pass(j) = 

end
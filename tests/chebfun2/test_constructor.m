function pass = test_constructor( pref ) 
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb2Prefs.chebfun2eps;

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
f5 = chebfun2( '2' );
f6 = chebfun2( @(x,y) 2 );
pass(4) = ( norm( f5 - f6 ) < tol && norm( f5 - 2 ) == 0 );

g = @(x,y) cos(x).*y + x.*sin(y); 
f = chebfun2( g );
pass(5) = ( length(f) == 2 ); 
pass(6) = ( abs( feval(f, .1, pi/6) - feval(g, .1, pi/6) ) < tol );

% George's function that failed:
f = @(x,y) 1./(1+25*x.^2.*y.^2);
ffch = chebfun2(@(x,y) f(x,y), [-2 2 -2 2]);
xx = linspace(-2,2); 
[XX,YY] = meshgrid(xx,xx);
pass(7) = ( max(max( abs(f(XX,YY) - ffch(XX,YY) ))) < 2e4*tol );

% Test that #928 is fixed.
f = chebfun2(@(x,y) 1, 'vectorize');
g = chebfun2(1);
pass(8) = norm(f - g) < tol;
f = chebfun2(@(x,y) 1, 'vectorize');
g = chebfun2(1);
pass(9) = norm(f - g) < tol;

% Test that the 2nd argument can be a preference:
p = pref;
p.tech = @trigtech;
f = chebfun2(@(x,y) sin(pi*x).*cos(pi*y), p);
g = chebfun2(@(x,y) sin(pi*x).*cos(pi*y), [-1 1 -1 1], p);
pass(10) = norm(f - g) < tol;

% Test the length of chebfun2s in different directions:
p.tech = @chebtech2;
f = chebfun2(@(x,y) sin(80*x+y), p);
pass(11) = length(f.cols(:,1)) < 50;

% Test whether there is a huge overesimation in the length of chebfun2s:
f = chebfun(@(x) 1./(1+25*x.^2));
f2 = chebfun2(@(x,y) 1./(1+25*x.^2));
pass(12) = length(f2.rows) < length(f)+20;

% Test whether reconstructing rows gives the same length.
f = chebfun2(@(x,y) cos(pi*x.*y));
frows = chebfun(f.rows);
pass(13) = abs(length(frows)-length(f.rows)) < 10;

end

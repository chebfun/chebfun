function pass = test_equiOption( pref ) 
% Test funqi in 2D 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb2Prefs.chebfun2eps;

% Canonical domain: 
dom = [-1 1 -1 1];
f = @(x,y) cos(x+y);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
[xx, yy] = meshgrid( x, y ); 
A = f(xx, yy) ; 
g = chebfun2( A , dom, 'equi' );
h = chebfun2( f, dom); 
pass(1) = norm( h - g ) < tol ; 

% Rectangular domain: 
dom = [-1 2 -2 1];
f = @(x,y) cos(x+y);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
[xx, yy] = meshgrid( x, y ); 
A = f(xx, yy) ; 
g = chebfun2( A , dom, 'equi' );
h = chebfun2( f, dom); 
pass(2) = norm( h - g ) < tol ; 

% Nonsymmetric function: 
dom = [-1 2 -2 1];
f = @(x,y) cos(x+2*y);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
[xx, yy] = meshgrid( x, y ); 
A = f(xx, yy) ; 
g = chebfun2( A , dom, 'equi' );
h = chebfun2( f, dom); 
pass(3) = norm( h - g ) < tol ; 

% Small domain; 
h = 1e-3;
dom = [1-h 1+h 1-2*h 1+2*h];
f = @(x,y) cos(x+2*y);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
[xx, yy] = meshgrid( x, y ); 
A = f(xx, yy) ; 
g = chebfun2( A , dom, 'equi' );
h = chebfun2( f, dom); 
pass(4) = norm( h - g ) < tol ; 

end 
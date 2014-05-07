function pass = test_basic_constructor( pref ) 
% This tests the chebfun2 constructor.  

if ( nargin < 1 ) 
    pref = chebfun2pref; 
end 
tol = 1e2 * pref.eps; 
j = 1; 

% Can we make a chebfun2: 
f = @(x,y) cos( x ) + sin( x .* y );  % simple function. 
fstr = 'cos(x) + sin(x.*y)'; % string version.

% construct with a domain.
f1 = chebfun2( f ); 
f2 = chebfun2( f, [-1 1 -1 1] );
pass(j) = ( norm(f1 - f2) < tol ); j = j+1; 

% construct from a string. 
f3 = chebfun2( fstr ); 
f4 = chebfun2( fstr, [-1 1 -1 1]); 
pass(j) = ( norm( f1 - f3 ) < tol ); j = j+1;
pass(j) = ( norm( f3 - f4 ) < tol ); j = j+1;

g = @(x,y) cos(x).*y + x.*sin(y); 
f = chebfun2( g );
pass(j) = ( length(f) == 2 ); j = j + 1; 
pass(j) = ( abs( feval(f, .1, pi/6) - feval(g, .1, pi/6) ) < tol ); j = j+1;


%%
% Constructing from matrix of values 
% A = randn( 100 ); 
% f = chebfun2( A ); 
% pass(j) = 

end
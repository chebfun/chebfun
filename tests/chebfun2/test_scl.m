function pass = test_scl( pref )
% Check correct vertical scaling. 


if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 10*pref.eps; 
j = 1;

% Scale invariant
f = chebfun2( @(x,y) cos(x.*y) );
g = chebfun2( @(x,y) eps*cos(x.*y) );
err = norm( eps*f - g);
pass(j) = ( err  <  tol ); j = j + 1; 

% hscale invariant 
f = chebfun2( @(x,y) cos(x.*y) );
g = chebfun2( @(x,y) cos(x/eps.*y/eps), eps*[-1 1 -1 1] );
err = abs( f(1,1) - g(eps,eps) );
pass(j) = ( err  < tol ); j = j + 1; 
err = abs( f(pi/6,1) - g(eps*pi/6,eps) );
pass(j) = ( err  < tol ); j = j + 1; 
end
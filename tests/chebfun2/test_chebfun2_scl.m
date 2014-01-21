function pass = test_chebfun2_scl( pref )
% Check correct vertical scaling. 


if ( nargin < 1 ) 
    pref = chebpref; 
end 
tol = pref.cheb2Prefs.eps; 
j = 1;

% Scale invariant
f = chebfun2( @(x,y) cos(x.*y) );
g = chebfun2( @(x,y) eps*cos(x.*y) );
pass(j) = ( norm( eps*f - g)  <  tol ); j = j + 1; 

% hscale invariant 
f = chebfun2( @(x,y) cos(x.*y) );
g = chebfun2( @(x,y) cos(x/eps.*y/eps), eps*[-1 1 -1 1] );
pass(j) = ( abs( f(1,1) - g(eps,eps) )  < eps* tol ); j = j + 1; 
pass(j) = ( abs( f(pi/6,1) - g(eps*pi/6,eps) )  < eps* tol ); j = j + 1; 
end
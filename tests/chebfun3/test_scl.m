function pass = test_scl( pref )
% Check correct vertical scaling. 


if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 10*pref.cheb3Prefs.chebfun3eps;
j = 1;

% Scale invariant
f = chebfun3(@(x,y,z) cos(x.*y.*z));
g = chebfun3(@(x,y,z) eps*cos(x.*y.*z));
err = norm( eps*f - g);
pass(j) = err  <  tol;
j = j + 1; 

% hscale invariant 
f = chebfun3(@(x,y,z) cos(x.*y.*z) );
g = chebfun3(@(x,y,z) cos(x/eps.*y/eps.*z/eps), eps*[-1 1 -1 1 -1 1]);
err = abs(f(1, 1, 1) - g(eps, eps, eps));
pass(j) = err < tol; 
j = j + 1;

err = abs(f(pi/6, 1, 1) - g(eps*pi/6, eps, eps) );
pass(j) = err  < tol;
j = j + 1; 
end
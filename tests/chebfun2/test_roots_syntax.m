function pass = test_roots_syntax( pref ) 
% Check the syntax to chebfun2/roots.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 10*pref.cheb2Prefs.chebfun2eps;

% Try out the many different ways of calling the same algorithms. 
f = chebfun2(@(x,y) cos(7*x.^2.*y + y)); 
g = chebfun2(@(x,y) cos(7*x.*y)); 
r1 = roots(f,g,'resultant');
r2 = roots(f,g,'ms'); 
r3 = roots(f,g,'marchingsquares');
r4 = roots(f,g);
r5 = roots([f;g],'resultant'); 
r6 = roots([f;g],'marchingsquares');
r7 = roots([f;g]);

pass(1) = ( norm( r1 - r5 ) == 0 );
pass(2) = ( norm( r2 - r3 ) == 0 );
pass(3) = ( norm( r2 - r6 ) == 0 );
pass(4) = ( norm( sort(r4) - sort(r6) ) < tol );
pass(5) = ( norm( sort(r4) - sort(r5) ) < tol );
pass(6) = ( norm( r4 - r7 ) == 0 );

end

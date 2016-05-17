function pass = test_cdr( ) 
% Test spherefun cdr() command

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% Pick a symmetric function:  
f = spherefun( @(lambda, theta) sin( theta ).*sin(lambda) );
[C, D, R] = cdr( f ); 
pass(1) = norm( C - R ) < tol;

% Check it works: 
f = spherefun( @(x,y,z) exp(-((x-.4).^2+(y-.9).^2+(z-.1).^2)) ); 
[C, D, R] = cdr( f );
pass(2) = norm( chebfun2(@(lam,th) feval(f,lam,th),[-pi,pi,-pi,pi]) - C * D * R' ) < tol;

end
function pass = test_cdr( ) 
% Test diskfun cdr() command

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;


% Check it works: 
f = diskfun( @(x,y,z) exp(-((x-.4).^2+(y-.9).^2))); 
[C, D, R] = cdr( f );
pass(1) = norm( chebfun2(@(t,r) feval(f,t,r, 'polar'),[-pi,pi,-1,1])...
- C * D * R' ) < tol;

end
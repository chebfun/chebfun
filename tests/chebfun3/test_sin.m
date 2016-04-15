function pass = test_sin( pref ) 
% This tests the sine function in Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

% This just constructs a function involving a sine:
f = chebfun3( @(x,y,z) sin(x+y+z) );
pass(1) = abs( f(.1,.2,.3) - sin(.6) ) < tol;

% This tests the Chebfun3 sine function:
f2 = sin(f);
pass(2) = abs( f2(.2,.2,.2) - sin(sin(.6)) ) < tol;

% This varies the domain:
d = [1 3 1 3 1 3];
x = chebfun3(@(x,y,z) x, d);
y = chebfun3(@(x,y,z) y, d);
z = chebfun3(@(x,y,z) z, d);
f3 = sin(2*x+3*y);
pass(3) = abs( f3(1,2,3) - sin(8) ) < tol;

% Here we check something in trig mode
f4 = chebfun3( @(x,y,z) sin(2*x + 3*y), [-pi pi -pi pi -pi pi], 'trig');
pass(4) = abs( sin(f4(1,1,1)) - sin(sin(5)) ) < tol;

ep = 1e-8;
tol2 = 1e2 * ep;
f5 = chebfun3( @(x,y,z) exp(x.*y.*z) , 'eps', ep);
pass(5) = abs( sin(f5(.5,.5,.5)) - sin(exp(0.125)) ) < tol2;

end

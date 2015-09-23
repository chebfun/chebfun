function pass = test_interpaccuracy( pref ) 
% Test the Chebfun2 constructor with a few functions.
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.eps; 

% Mohsin's bugs. 
f = @(x,y) exp(pi*(x+y)); 
x=.995; y=.99; 
g = chebfun2(f); 
pass(1) = (abs(f(x,y)-g(x,y))<3000*tol);  

% composing the same function. 
x = chebfun2(@(x,y)x); y = chebfun2(@(x,y)y); 
f = exp(pi*(x+y)); 
x=.995; y=.99;
pass(2) = (abs(f(x,y)-g(x,y))<1e4*tol);


% Nick T's function 
f = @(x,y) exp(-100*( x.^2 - x.*y + 2*y.^2 - 1/2).^2);
g = chebfun2(f);
twonorm = 0.545563608722019;
pass(3) = (abs(norm(g)-twonorm)<100*tol); 

x=.995; y=.99;
pass(4) = (abs(f(x,y)-g(x,y))<100*tol);  


f = chebfun2(@(x,y) exp(3*(cos(x-.1)+cos(y-.2))));
pass(5) = (abs(f(.1,.2) - exp(6))<1e4*tol); 

end
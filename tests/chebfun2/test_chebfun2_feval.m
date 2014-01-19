function pass = test_chebfun2_feval( pref ) 
% Evaluation check for Chebfun2. 

if ( nargin < 1 ) 
    pref = chebpref; 
end 

tol = 100*pref.cheb2Prefs.eps; 
j = 1; 

f = chebfun2(@(x,y) x); 
pass(j) = (abs(f(pi/6,pi/12)-pi/6) < tol); j=j+1; 

f = chebfun2(@(x,y) x, [-1 2 -pi/2 pi]); 
pass(j) = (abs(f(0,0)) < 1e-14);  j=j+1; 
pass(j) = (abs(f(pi/6,pi/12)-pi/6) < tol);  j=j+1; 


f = chebfun2(@(x,y) y, [-1 2 -pi/2 pi]); 
pass(j) = (abs(f(0,0)) < 1e-14);  j=j+1; 
pass(j) = (abs(f(pi/6,pi/12)-pi/12) < tol);  j=j+1; 


% some harder tests. 
f = @(x,y) cos(x) + sin(x.*y); 
g = chebfun2(f);

r = 0.126986816293506; s = 0.632359246225410; % two fixed random number in domain.
pass(j) = (abs(f(r,s) - g(r,s))<tol);j=j+1;

% Are we evaluating on arrays correctly
r = rand(10,1); s = rand(10,1); [rr, ss]=meshgrid(r,s);
pass(j) = (norm((f(r,s) - g(r,s)))<tol);j=j+1;
pass(j) = (norm((f(rr,ss) - g(rr,ss)))<tol); j=j+1;% on arrays as well. 

% Does this work off [-1,1]^2
g = chebfun2(f,[-pi/6 pi/2 -pi/12 sqrt(3)]); % strange domain. 
r = 0.126986816293506; s = 0.632359246225410; % two fixed random number in domain.
pass(j) = (abs(f(r,s) - g(r,s))<tol);j=j+1;

% Are we evaluating on arrays correctly
r = rand(10,1); s = rand(10,1); [rr ss]=meshgrid(r,s);
pass(j) = (norm((f(r,s) - g(r,s)))<tol);j=j+1;
pass(j) = (norm((f(rr,ss) - g(rr,ss)))<tol);j=j+1; % on arrays as well. 

% Evaluation at complex arguments. 
f = chebfun2(@(x,y) x+1i*y); 
pass(j) = ( norm( f(-1,-1) - (-1-1i) )  < tol ); j = j + 1; 
pass(j) = ( norm( f(r,s) - (r + 1i*s) )  < tol );j = j + 1; 
pass(j) = (norm((f(rr,ss) - (rr + 1i*ss)))<tol); j=j+1; % on arrays as well. 

% Evaluation at complex arguments, different syntax. 
f = chebfun2(@(z) z); 
pass(j) = ( norm( f(-1-1i) - (-1-1i) )  < tol ); j = j + 1; 
pass(j) = ( norm( f(r+1i*s) - (r + 1i*s) )  < tol );j = j + 1; 
pass(j) = (norm((feval(f,rr,ss) - (rr + 1i*ss)))<tol); j=j+1; % on arrays as well. 


end
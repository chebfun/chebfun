function pass = test_feval( pref ) 
% Test feval

if ( nargin == 0) 
    pref = chebfunpref; 
end

seedRNG(42);

tol = 100*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x, [-1 2 -pi/2 pi -3 1]); 
pass(1) = (abs(f(0,0,0)) < tol*vscale(f));  
pass(2) = (abs(f(pi/6,pi/12,-1)-pi/6) < tol*vscale(f));  

f = chebfun3(@(x,y,z) y, [-1 2 -pi/2 pi -3 1]); 
pass(3) = (abs(f(0,0,0)) < tol);   
pass(4) = (abs(f(pi/6,pi/12,-1)-pi/12) < tol*vscale(f)); 

f = chebfun3(@(x,y,z) z, [-1 2 -pi/2 pi -3 1]); 
pass(5) = (abs(f(0,0,0)) < tol);   
pass(6) = (abs(f(pi/6,pi/12,-1)+1) < tol*vscale(f)); 

% some harder tests. 
f = @(x,y,z) cos(x) + sin(x.*y) + sin(z.*x); 
g = chebfun3(f);

pts = 2*rand(3,1) - 1;
pass(7) = (abs(f(pts(1),pts(2),pts(3)) - g(pts(1),pts(2),pts(3)))<tol*vscale(g));

% Are we evaluating on arrays correctly
r = rand(10,1); 
s = rand(10,1); 
t = rand(10,1); 
[rr, ss, tt]=meshgrid(r,s,t);
pass(8) = (norm((f(r,s,t) - g(r,s,t))) < tol*vscale(g));
pass(9) = (max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt))))) < tol*vscale(g));

% Does this work off [-1,1]^2
g = chebfun3(f,[-pi/6 pi/2 -pi/12 sqrt(3) -3 1]); % strange domain. 
r = 0.126986816293506; s = 0.632359246225410; t = 0.351283361405006;
% three fixed random number in domain.
pass(10) = (abs(f(r,s,t) - g(r,s,t))<tol*vscale(g));

% Are we evaluating on arrays correctly
pass(11) = (norm((f(r,s,t) - g(r,s,t)))<tol*vscale(g));
pass(12) = (max(max(max(abs(f(rr,ss,tt) - g(rr,ss,tt)))))<tol*vscale(g)); 

% [TODO] trig, meshgrid, ndgrid?

end

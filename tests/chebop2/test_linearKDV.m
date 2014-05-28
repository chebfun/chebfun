function pass = chebop2_linearKDV
% Check that we can solve linear KDV equations. 
% Alex Townsend, April 2013. 

tol = 100*(chebfun2pref('eps'));
j = 1; 

% Simple example. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) diff(u)-exp(-t).*exp(1)];
N.lbc = @(t) exp(-t).*exp(-1);
u = N \ 0;
 
pass(j) = ( norm(u-exact) < tol); j = j + 1; 


% Another simple example. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) diff(u)-exp(-t).*exp(1)];
N.lbc = @(t,u) diff(u) - exp(-t).*exp(-1);
u = N \ 0;
 
pass(j) = ( norm(u-exact) < tol); j = j + 1; 

% Different boundary conditions. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) diff(u,2)-exp(-t).*exp(1)];
N.lbc = @(t,u) diff(u) - exp(-t).*exp(-1);
u = N \ 0;
 
pass(j) = ( norm(u-exact) < tol); j = j + 1; 

% Different boundary conditions. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) diff(u,2)-exp(-t).*exp(1)];
N.lbc = @(t,u) u - exp(-t).*exp(-1);
u = N \ 0;
 
pass(j) = ( norm(u-exact) < tol); j = j + 1; 


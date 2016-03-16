function pass = test_linearKDV( pref )
% Check that we can solve linear KDV equations. 
% Alex Townsend, April 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 100*pref.cheb2Prefs.chebfun2eps;

% Simple example. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) ; diff(u)-exp(-t).*exp(1)];
N.lbc = @(t) exp(-t).*exp(-1);
u = N \ 0;
 
pass(1) = ( norm(u-exact) < tol); 


% Another simple example. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) ; diff(u)-exp(-t).*exp(1)];
N.lbc = @(t,u) diff(u) - exp(-t).*exp(-1);
u = N \ 0;
 
pass(2) = ( norm(u-exact) < 2*tol); 

% Different boundary conditions. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) ; diff(u,2)-exp(-t).*exp(1)];
N.lbc = @(t,u) diff(u) - exp(-t).*exp(-1);
u = N \ 0;
 
pass(3) = ( norm(u-exact) < 300*tol);

% Different boundary conditions. 
d = [-1 1 0 1];
exact = chebfun2(@(x,t) exp(-t).*exp(x),d);

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(x);
N.rbc = @(t,u) [u - exp(-t).*exp(1) ; diff(u,2)-exp(-t).*exp(1)];
N.lbc = @(t,u) u - exp(-t).*exp(-1);
u = N \ 0;
 
pass(4) = ( norm(u-exact) < 100*tol); 


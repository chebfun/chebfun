function pass = test_squarewaveequation( pref )
% Testing wave equation on a square domain. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 100*pref.cheb2Prefs.chebfun2eps;

%%
% First example. 
d = [0 10 0 10]; 
exact = chebfun2(@(x,t) sin(x+t) + sin(x-t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(t) + sin(-t);
N.rbc = @(t) sin(10+t) + sin(10-t);
N.dbc = @(x,u) [u - 2*sin(x) ; diff(u)];
u = N \ 0;

pass(1) = ( norm(u - exact) < 150*tol ); 

%%
d = [0 1 0 1]; 
exact = chebfun2(@(x,t) sin(x+t) + cos(x-t),d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2),d);
N.lbc = @(t) sin(t) + cos(-t);
N.rbc = @(t) sin(1+t) + cos(1-t);
N.dbc = @(x,u) [u - sin(x) - cos(x) ; diff(u)-cos(x)-sin(x)];
u = N \ 0;

pass(2) = ( norm(u - exact) < tol ); 

%%
d = [0 pi 0 pi]; 
exact = chebfun2(@(x,t) sin(x+t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(t);
N.rbc = @(t) sin(pi+t);
N.dbc = @(x,u) [u - sin(x) ; diff(u) - cos(x)];
u = N \ 0;

pass(3) = ( norm(u - exact) < 5*tol ); 


%% different wave number. 
% d = [0 pi 0 pi]; 
% exact = chebfun2(@(x,t) sin(x+2*t), d); 
% N = chebop2(@(u) diff(u,2,1) - 2*diff(u,2,2), d);
% N.lbc = @(t) sin(2*t);
% N.rbc = @(t) sin(pi+2*t);
% N.dbc = @(x,u) [u - sin(x) diff(u) - 2*cos(x)];
% u = N \ 0;
% 
% pass(j) = ( norm(u - exact) < tol ); j = j + 1 ;



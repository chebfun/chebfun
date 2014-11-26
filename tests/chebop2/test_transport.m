function pass = test_transport( prefs )
% Testing the transport equation on rectangular domains.
% Alex Townsend, April 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.eps; 

%% Simple example of a transport equation. 

d = [-1,1,0,1];
exact = chebfun2(@(x,t) exp(-x).*exp(t),d);
N = chebop2(@(u) diffy(u,1) + diffx(u,1), d);
N.dbc = @(x) exp(-x);
N.lbc = @(t) exp(1).*exp(t);
u = N \ 0;
pass(1) = ( abs(norm(u - exact) < 2*tol)); 

%% Simple example on square domain. 
d = [0,1,0,1];
exact = chebfun2(@(x,t) exp(-x).*exp(t),d);
N = chebop2(@(u) diffy(u,1) + diffx(u,1), d);
N.dbc = @(x) exp(-x);
N.lbc = @(t) exp(t);
u = N \ 0;
pass(2) = ( abs(norm(u - exact) < tol)); 

%% Simple example on square domain. 
d = [-1,1,0,1];
exact = chebfun2(@(x,t) exp(-x).*exp(t) + exp(-.5*x).*exp(.5*t), d);
N = chebop2(@(u) diffy(u,1) + diffx(u,1), d);
N.dbc = @(x) exp(-x) + exp(-.5*x);
N.lbc = @(t) exp(1).*exp(t) + exp(.5).*exp(.5*t);
u = N \ 0;
pass(3) = ( abs(norm(u - exact) < 5*tol)); 


%% Transport equation with different transport parameter. 
d = [-pi,pi,0,1];
exact = chebfun2(@(x,t) exp(x).*exp(-5*t), d);
N = chebop2(@(u) diffy(u,1) + 5*diffx(u,1), d);
N.dbc = @(x) exp(x);
N.lbc = @(t) exp(-pi).*exp(-5*t);
u = N \ 0;
pass(4) = ( norm(u - exact) < 10*tol );  

%% Transport equation with different transport parameter, and large time
% length
d = [-pi,pi,0,100];
exact = chebfun2(@(x,t) exp(x).*exp(-t/10), d);
N = chebop2(@(u) diffy(u,1) + .1*diffx(u,1), d);
N.dbc = @(x) exp(x);
N.lbc = @(t) exp(-pi).*exp(-t/10);
u = N \ 0;
pass(5) = ( norm(u - exact) < 100*tol);


% %% Here's another one. 
% % Just check that the code runs, without checking error. 
% d = [0,1,0,1];
% N = chebop2(@(u) diffy(u,1) + diffx(u,1), d);
% N.dbc = @(x) sech(10*x);
% N.lbc = @(t) sech(10*t);
% u = N \ 0;
% x = chebfun('x',[0,1]);
% pass(j) = ( abs(norm(u(:,0) - sech(10*x)) < tol)); j = j + 1;
% pass(j) = ( abs(norm(u(0,:) - sech(10*x)) < tol)); j = j + 1;




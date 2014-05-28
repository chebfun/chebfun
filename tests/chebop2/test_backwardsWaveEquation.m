function pass = chebop2_backwardsWaveEquation
% Check that the backwards wave equation is working. 
% Alex Townsend, August 2013. 

j = 1; 
tol = chebfun2pref('eps');

%%

d = [-pi pi 0 1]; 
exact = chebfun2(@(x,t) sin(x+t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(-pi+t);
N.rbc = @(t) sin(pi+t);
N.ubc = @(x,u) [u - sin(x+1) diff(u) - cos(x+1)];
u = N \ 0;
pass(j) = ( norm(u - exact) < tol ) ; j = j + 1; 


%%
% d = [-pi pi 0 1]; 
% exact = chebfun2(@(x,t) exp(x+1i*t), d); 
% N = chebop2(@(u) diff(u,2,1) + diff(u,2,2), d);
% N.lbc = @(t) exp(-pi+1i*t);
% N.rbc = @(t) exp(pi+1i*t);
% N.ubc = @(x,u) [u - exp(x+1i) diff(u) - 1i*exp(x+1i)];
% u = N \ 0;
% x = chebpts(100,d(1:2)); y = chebpts(100,d(3:4));
% [xx,yy]=meshgrid(x,y);
% 
% pass(j) = (norm(u(xx,yy) - exact(xx,yy),inf) <5e-3); j = j +1; 

end
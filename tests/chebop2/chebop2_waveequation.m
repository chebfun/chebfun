function pass = chebop2_waveequation
% Testing the wave equation on a non-square domain.  This tests if the domain is
% being treated properly and if Neumann boundary conditions are being dealt with
% correctly. 
% Alex Townsend, April 2013. 

tol = 100*chebfun2pref('eps');
j = 1; 

%%  Standard wave equation on a non-square domain. 

d = [-pi pi 0 1]; 
exact = chebfun2(@(x,t) sin(x+t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(-pi+t);
N.rbc = @(t) sin(pi+t);
N.dbc = @(x,u) [u - sin(x) diff(u) - cos(x)];
u = N \ 0;
pass(j) = (abs( norm(u - exact)) < tol); j = j + 1; 

%% Another standard example on a non-square domain. 

d = [-pi pi 0 1]; 
exact = chebfun2(@(x,t) sin(x+t) + sin(x-t), d); 
N = chebop2(@(u) diff(u,2,1) - diff(u,2,2), d);
N.lbc = @(t) sin(-pi+t) + sin(-pi -t);
N.rbc = @(t) sin(pi+t) + sin(pi -t);
N.dbc = @(x,u) [u - 2*sin(x) diff(u)];
u = N \ 0;
pass(j) = (abs( norm(u - exact)) < tol); j = j + 1;

%% An example with a different wave speed. 

d = [-2*pi 2*pi 0 1]; c = 2; 
exact = chebfun2(@(x,t) sin(x+c*t) + sin(x-c*t), d); 
N = chebop2(@(u) diff(u,2,1) - c.^2*diff(u,2,2), d);
N.lbc = @(t) sin(-2*pi+c*t) + sin(-2*pi - c*t);
N.rbc = @(t) sin(2*pi+c*t) + sin(2*pi - c*t);
N.dbc = @(x,u) [u - 2*sin(x) diff(u)];
u = N \ 0;
pass(j) = (abs( norm(u - exact)) < tol); j = j + 1;

%% Higher wave speed. 

d = [-2*pi 2*pi 0 1]; c = 30; 
exact = chebfun2(@(x,t) sin(x+c*t), d); 
N = chebop2(@(u) diff(u,2,1) - c^2*diff(u,2,2), d);
N.lbc = @(t) sin(-2*pi+c*t);
N.rbc = @(t) sin(2*pi+c*t);
N.dbc = @(x,u) [u - sin(x) diff(u) - c*cos(x)];
u = N \ 0;
[xx,yy] = meshgrid(linspace(d(1),d(2)),linspace(d(3),d(4)));
pass(j) = ( norm(u(xx,yy) - exact(xx,yy), inf) < 1e5*tol); j = j + 1;

%% Working for non-zero starting time. 
d = [-2*pi 2*pi 1 2]; c = 3; 
exact = chebfun2(@(x,t) sin(x+c*t), d); 
N = chebop2(@(u) diff(u,2,1) - c^2*diff(u,2,2), d);
N.lbc = @(t) sin(-2*pi+c*t);
N.rbc = @(t) sin(2*pi+c*t);
N.dbc = @(x,u) [u - sin(c+x) diff(u) - c*cos(x+c)];
u = N \ 0;
pass(j) = (abs( norm(u - exact)) < tol); j = j + 1;

end
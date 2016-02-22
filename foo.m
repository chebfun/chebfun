%dom = [0,10];
%a = 1;
%b = 1;
%alpha = 1;
%beta = 0;
%x = chebfun('x', dom);

%G = chebop(@(x,u,w) .5*(a*u.^2 + b*w.^2), dom);

%F = chebop(@(x,u,w) diff(u,2) + sin(u) - w, dom);
%F.lbc = @(u,w) [u-alpha; diff(u)-beta]; 

%u0 = alpha + beta*(x-1);
%u0 = 0*x;
%w0 = 0*x;

%[u,w,Gval] = fmincon(G,F,u0,w0);



dom = [-1,1];
a = 1;
b = 10;
alpha = 1;
beta = 1.3;
x = chebfun('x', dom);

G = chebop(@(x,u,w) .5*(a*u.^2 + b*w.^2), dom);

F = chebop(@(x,u,w) .2*diff(u,3) + cos(u) - w, dom);
F.lbc = @(u,w) [u-alpha;diff(u)]; 
F.rbc = @(u,w) u-beta; 

u0 = 0*x;
w0 = 0*x;

[u,w,Gval] = fmincon(G,F,u0,w0);



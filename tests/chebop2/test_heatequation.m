function pass = test_heatequation( prefs ) 
% Check that the solver can do heat equation. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 1e10*prefs.techPrefs.eps; 

% Solution to heat equation should not depend on the time interval. 
d = [-1 1 0 1]; k = 1; 
N = chebop2(@(u) diff(u,1,1) - k*diff(u,2,2), d);
N.dbc = @(x) -(x+1).*(x-1); N.lbc = 0; N.rbc = 0; 
u = N \ 0;

d = [-1 1 0 1.5]; k = 1;  
N = chebop2(@(u) diff(u,1,1) - k*diff(u,2,2), d);
N.dbc = @(x) -(x+1).*(x-1); N.lbc = 0; N.rbc = 0; 
v = N \ 0;

d = [-1 1 0 1];
x = chebpts(100,d(1:2));
y = chebpts(100,d(3:4)); 
[xx, yy] = meshgrid(x,y); 

% res = v{d(1),d(2),d(3),d(4)}; 
pass(1) = ( norm( u(xx,yy) - v(xx,yy), inf) < 500*tol ); 

% Do we agree with pde15s on a simple example? 
k = 2;  
f = chebfun2(@(x,t) exp(-40*x.^2),[-1 1 0 1]);

% chebfun/pde15s
bc.left = 0;
bc.right = 0;
opts = pdeset('plot','off');
uu = pde15s(@(x,t,u) k*diff(u,2),0:.1:1,f(:,0),bc,opts);

% chebop2 
d = [-1 1 0 1];
N = chebop2(@(u) diff(u,1,1) - k*diff(u,2,2), d);
N.dbc = f(:,0); N.lbc = 0; N.rbc = 0; 
u = N \ 0;

pass(2) = ( norm(uu(:,end) - u(:,1).') < tol ); 

% Do we agree with pde15s on another simple example? 
% Different heat coefficient.
d = [-1 1 0 1];
k = pi;  
f = chebfun2(@(x,t) exp(-40*x.^2),d);

% chebfun/pde15s
bc.left = 0;
bc.right = 0;
opts = pdeset('plot','off');
uu = pde15s(@(x,t,u) k*diff(u,2),0:.1:1,f(:,0),bc,opts);

% chebop2 
N = chebop2(@(u) diff(u,1,1) - k*diff(u,2,2), d);
N.dbc = f(:,0); N.lbc = 0; N.rbc = 0; 
u = N \ 0;

pass(3) = ( norm(uu(:,end) - u(:,1).') < tol ); 

% Do we agree with pde15s on another simple example? 
% Different domain in space.
d = [-2 2 0 1];
k = 1;  
f = chebfun2(@(x,t) exp(-40*x.^2),d);

% chebfun/pde15s
bc.left = 0;
bc.right = 0;
opts = pdeset('plot','off');
uu = pde15s(@(x,t,u) k*diff(u,2),0:.1:1,f(:,0),bc,opts);

% chebop2 
N = chebop2(@(u) diff(u,1,1) - k*diff(u,2,2), d);
N.dbc = f(:,0); N.lbc = 0; N.rbc = 0; 
u = N \ 0;

pass(4) = ( norm(uu(:,end) - u(:,1).') < tol );

end

% % Heat equation 
% N = chebop2(@(u) diff(u,1,1) - .1*diff(u,2,2), [-4 4 0 6]); 
% N.dbc = chebfun(@(x) sin(pi*x/4) + sin(pi*3*x/4) + sin(pi*12*x/4), [-4,4]); 
% N.lbc = 0; 
% N.rbc = 0; 
% u = N \ 0; 
% 
% for t = [0:.001:.2]
%    plot(u(:,t)), hold on, shg
%    axis([-4 4 -3 3]), pause(.01)
%    hold off
% end

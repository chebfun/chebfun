function pass = test_advectionDiffusion1( prefs ) 
% Check advection diffusion equations. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 1e10*prefs.techPrefs.eps; 

% Advection-diffusion 1
N = chebop2(@(u) diff(u,1,1) - .1*diff(u,2,2) - diff(u,1,2), [-2.5 3 0 6]); 
N.dbc = chebfun(@(x) sin(pi*x), [-2.5 3]); 
N.lbc = @(x,u) diff(u);  
N.rbc = 0; 
u = N \ 0;

% chebfun/pde15s
f = chebfun(@(x) sin(pi*x), [-2.5 3]);
bc.left = @(u) diff(u);
bc.right = 0;
opts = pdeset('plot','off');
uu = pde15s(@(t,x,u) .1*diff(u,2) + diff(u), 0:.1:6, f, bc, opts);

k = 1; j = 1; 
for t = 0:.1:6
   pass(j) = ( norm(u(:,t) - uu(:,k).') < 10*tol ); 
   j = j + 1; 
   k = k + 1; 
end

end
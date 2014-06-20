function pass = test_advectionDiffusion2( prefs ) 
% Check advection diffusion equations. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 1e10*prefs.techPrefs.eps; 

% Advection-diffusion 2
N = chebop2(@(u) diff(u,1,1) - .3*diff(u,2,2) - 10*diff(u,1,2), [-1 1 0 .25]); 
N.dbc = chebfun(@(x) exp(-10*x.^4./(1-x.^2))); 
N.lbc = 0; 
N.rbc = 0; 
u = N \ 0; 

% chebfun/pde15s
f = chebfun(@(x) exp(-10*x.^4./(1-x.^2)),[-1 1]);
bc.left = 0;
bc.right = 0; 
opts = pdeset('plot','off');
uu = pde15s(@(t,x,u) .3*diff(u,2) + 10*diff(u), 0:.005:.25, f, bc, opts);

k = 1; j = 1;
for t = 0:.005:.25
   pass(j) = ( norm(u(:,t) - uu(:,k).') < 2*tol ); 
   j = j + 1;
   k = k + 1; 
end

end

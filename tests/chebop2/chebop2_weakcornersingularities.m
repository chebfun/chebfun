function pass = chebop2_weakcornersingularities
% Check that we can correctly resolve weak corner singularities to high
% accuracy. 
% Alex Townsend, August 2013. 

tol = 1e6*chebfun2pref('eps'); 
j = 1; 

f = @(x,y) (x.^2+y.^2).*(log(sqrt(x.^2+y.^2)).*sin(2*atan(y./x)) + atan(y./x).*cos(2*atan(y./x)));

e = 1e-16;   % take us away from the corner a little (o/w f is NaN). 
N = chebop2(@(u) diffx(u,2) + diffy(u,2), [e 1-e e 1-e]); 
N.lbc = chebfun(@(y) f(e,y), [e 1-e]); 
N.rbc = chebfun(@(y) f(1-e,y), [e 1-e]); 
N.dbc = chebfun(@(x) f(x,e), [e 1-e]); 
N.ubc = chebfun(@(x) f(x,1-e), [e 1-e]); 
u = N \ 0; 
g = chebfun2(f, [e 1-e e 1-e]);

x = pi/6; y = pi/12; 
pass(j) = ( abs(g(x,y) - u(x,y)) < tol ); j = j +1; 
x = .8; y = .23; 
pass(j) = ( abs(g(x,y) - u(x,y)) < tol ); j = j +1; 

end

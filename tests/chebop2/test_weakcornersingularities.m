function pass = test_weakcornersingularities( pref )
% Check that we can correctly resolve weak corner singularities to high
% accuracy. 
% Alex Townsend, August 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 1e6*pref.cheb2Prefs.chebfun2eps;

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
pass(1) = ( abs(g(x,y) - u(x,y)) < tol );
x = .8; y = .23; 
pass(2) = ( abs(g(x,y) - u(x,y)) < tol ); 

end

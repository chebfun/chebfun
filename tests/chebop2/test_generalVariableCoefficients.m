function pass = test_generalVariableCoefficients( prefs )
% Check that we are solving general variable coefficient PDEs correctly. 
% Alex Townsend, August 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 4000*prefs.cheb2Prefs.chebfun2eps;

x = chebfun2(@(x,y) x); y = chebfun2(@(x,y) y); 
m = cos(x.*y);
exact = chebfun2(@(x,y) sin(5*x.*y)); 
N = chebop2(@(x,y,u) diff(u,2,2) - m.*diff(u,2,1));
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.dbc = exact(:,-1); N.ubc = exact(:,1); 
u = N \ (diff(exact,2,2) - m.*diff(exact,2,1));
err = norm(u - exact);
pass(1) = ( err < tol );


m = exp(-x).*y + y.^2;
exact = chebfun2(@(x,y) sin(5*x.*y)); 
N = chebop2(@(x,y,u) diff(u,2,2) - m.*diff(u,2,1));
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.dbc = exact(:,-1); N.ubc = exact(:,1); 
u = N \ (diff(exact,2,2) - m.*diff(exact,2,1));
err = norm(u - exact);
pass(2) = ( err < tol ); 


m1 = exp(-x).*y + y.^2; m2 = cos(x) + sin(y) + x.*y;
exact = chebfun2(@(x,y) sin(x.*y)); 
N = chebop2(@(x,y,u) m1.*diff(u,2,2) - m2.*diff(u,2,1));
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.dbc = exact(:,-1); N.ubc = exact(:,1); 
u = N \ (m1.*diff(exact,2,2) - m2.*diff(exact,2,1));
err = norm(u - exact);
pass(3) = ( err < tol ); 


m1 = exp(-x).*y + y.^2; m2 = cos(x) + sin(y) + x.*y;
exact = chebfun2(@(x,y) exp(-x).*y + sin(x) + cos(x.*y)); 
N = chebop2(@(x,y,u) m1.*diff(u,2,2) - m2.*diff(u,2,1));
N.lbc = exact(-1,:); N.rbc = exact(1,:); 
N.dbc = exact(:,-1); N.ubc = exact(:,1); 
u = N \ (m1.*diff(exact,2,2) - m2.*diff(exact,2,1));
err = norm(u - exact);
pass(4) = ( tol < 10*tol ); 

d = [-2 2 -2 2];
x = chebfun2(@(x,y) x, d); y = chebfun2(@(x,y) y, d); 
m = cos(x);
exact = chebfun2(@(x,y) sin(x.*y), d); 
N = chebop2(@(x,y,u) diff(u,2,2) - m.*diff(u,2,1), d);
N.lbc = exact(d(1),:); N.rbc = exact(d(2),:); 
N.dbc = exact(:,d(3)); N.ubc = exact(:,d(4)); 
u = N \ (diff(exact,2,2) - m.*diff(exact,2,1));
err = norm(u - exact);
pass(5) = ( err < 5*tol ); 

d = [-2 0 -4.1 pi];
x = chebfun2(@(x,y) x, d); y = chebfun2(@(x,y) y, d); 
m = x+y;
exact = chebfun2(@(x,y) 1+x+y+x.*y+x.^2.*y.^2, d); 
N = chebop2(@(x,y,u) diff(u,2,2) - m.*diff(u,2,1), d);
N.lbc = exact(d(1),:); N.rbc = exact(d(2),:); 
N.dbc = exact(:,d(3)); N.ubc = exact(:,d(4)); 
u = N \ (diff(exact,2,2) - m.*diff(exact,2,1));
err = norm(u - exact);
pass(6) = ( err < tol ); 


end
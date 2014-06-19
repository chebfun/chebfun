function pass = test_rhs( prefs )
% Check that non-zero righthand forcing terms are being dealt with
% correctly. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.eps; 

d = [0 pi 0 pi]; 
N = chebop2(@(u) laplacian(u), d); 
N.lbc = 0; 
N.rbc = @(y) pi*y.^3./6; 
N.dbc = 0; 
N.ubc = @(x) x*pi^3/6 + sin(x).*sinh(pi); 

rhs = chebfun2(@(x,y) x.*y, d); 

u = N \ rhs;

exact = chebfun2(@(x,y) x.*y.^3/6 + sin(x).*sinh(y), d);

pass(1) = ( norm(u - exact) < 10*tol );  



d = [0 pi 0 pi]; 
N = chebop2(@(u) laplacian(u), d); 
N.lbc = 0; 
N.rbc = @(y) pi^4*y./12 + sinh(pi)*(cos(y) + sin(y)); 
N.dbc = @(x) sinh(x); 
N.ubc = @(x) pi*x.^4/12 - sinh(x); 

rhs = chebfun2(@(x,y) x.^2.*y, d); 

u = N \ rhs;

exact = chebfun2(@(x,y) x.^4.*y./12 + sinh(x).*(cos(y) + sin(y)), d);

pass(2) = ( norm(u - exact) < 10*tol );


% right hand side when we are on a rectangular domain. 
% real part of analytic plus forcing term
d = [0 pi 0 1];
exact = chebfun2(@(x,y) real(exp(x+1i*y)) + x.^3.*y.^3, d);
N = chebop2(@(u) laplacian(u), d); 
N.lbc = exact(d(1),:); 
N.rbc = exact(d(2),:); 
N.dbc = exact(:,d(3)); 
N.ubc = exact(:,d(4)); 

u = N \ chebfun2(@(x,y) 6*x.*y.^3  + 6*y.*x.^3, d); 

pass(3) = ( norm(u - exact) < 10*tol );

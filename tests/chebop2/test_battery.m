function pass = test_battery( prefs ) 
% Check that chebop2 is working by using a battery of Laplace problems. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.eps; 

% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real(exp(x+1i*y)); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(1) = ( norm( exact - u ) < tol ); 

% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real(exp(2*(x+1i*y))); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(2) = ( norm( exact - u ) < tol ); 


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) 10*real(exp(2*(x+1i*y))); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(3) = ( norm( exact - u ) < 10*tol ); 


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real((x+1i*y).^2); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(4) = ( norm( exact - u ) < tol ); 


% Linearity check; 
M = N; 

N.lbc = @(x) (1+x).*(1-x); 
N.rbc = 0; N.ubc = 0; N.dbc = 0; 
u1 = N \ 0; 

N = M; N.rbc = @(x) (1+x).*(1-x); 
N.lbc = 0; N.ubc = 0; N.dbc = 0; 
u2 = N \ 0; 

N = M; N.ubc = @(x) (1+x).*(1-x); 
N.rbc = 0; N.lbc = 0; N.dbc = 0; 
u3 = N \ 0; 

N = M; N.dbc = @(x) (1+x).*(1-x); 
N.rbc = 0; N.ubc = 0; N.lbc = 0; 
u4 = N \ 0; 

N = M; 
N.lbc = @(x) (1+x).*(1-x);
N.rbc = @(x) (1+x).*(1-x);  
N.dbc = @(x) (1+x).*(1-x); 
N.ubc = @(x) (1+x).*(1-x);
u = N \ 0; 

[xx, yy] = chebfun2.chebpts2(100); 
A = feval(u,xx,yy); 
B = feval(u1,xx,yy) + feval(u2,xx,yy) + feval(u3,xx,yy) + feval(u4,xx,yy); 
pass(5) = ( norm( A - B ) < 1e10*tol);

% Check we can use the notation lap(u) = div(grad(u))
N = chebop2(@(u) -divergence(gradient(u)) );
N.lbc = 0; N.rbc = 0; N.dbc = 0; N.ubc = 0; 
u = N \ 1;
N = chebop2(@(u) -laplacian(u) );
N.lbc = 0; N.rbc = 0; N.dbc = 0; N.ubc = 0; 
exact = N \ 1;
pass(6) = ( norm( u - exact ) < tol);
end

function pass = chebop2_battery
% Check that chebop2 is working by using a battery of Laplace problems. 
% Alex Townsend, March 2013. 

j = 1; 
tol = chebfun2pref('eps');

% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real(exp(x+1i*y)); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(j) = ( norm( exact - u ) < tol ); j = j + 1; 


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real(exp(2*(x+1i*y))); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(j) = ( norm( exact - u ) < tol ); j = j + 1; 


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) 10*real(exp(2*(x+1i*y))); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(j) = ( norm( exact - u ) < tol ); j = j + 1; 


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
bdy = @(x,y) real((x+1i*y).^2); 
N.lbc = @(y) bdy(-1,y); N.rbc = @(y) bdy(1,y); 
N.dbc = @(x) bdy(x,-1); N.ubc = @(x) bdy(x,1); 
u = N \ 0; exact = chebfun2(bdy);
pass(j) = ( norm( exact - u ) < tol ); j = j + 1;


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

[xx yy] = chebpts2(100); 
A = feval(u,xx,yy); 
B = feval(u1,xx,yy) + feval(u2,xx,yy) + feval(u3,xx,yy) + feval(u4,xx,yy); 
pass(j) = ( norm( A - B ) < 1e10*tol); j = j + 1; 


end

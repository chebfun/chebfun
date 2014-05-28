function pass = chebop2_domain
%Check that boundary conditions are being imposed correctly. 
% Alex Townsend, March 2013. 

tol = 100*chebfun2pref('eps'); 
j = 1; 

d = [-2 2 -2 2];
N = chebop2(@(u) diff(u, 2, 1) + diff(u, 2, 2), d);

N.rbc = @(x) (2-x).*(2+x); 
N.ubc = 0; 
N.dbc = 0; 
N.lbc = 0;

u = N \ 0; 

pass(j) = ( norm( u(:,d(3)) - N.dbc ) < tol); j = j + 1; 
pass(j) = ( norm( u(:,d(4)) - N.ubc ) < tol); j = j + 1;
pass(j) = ( norm( u(d(1),:) - N.lbc ) < tol); j = j + 1;
pass(j) = ( norm( u(d(2),:) - N.rbc ) < tol); j = j + 1; 


d = [-pi 2*pi -2 5];
N = chebop2(@(u) diff(u, 2, 1) + diff(u, 2, 2), d);

N.rbc = 0; 
N.ubc = 0; 
N.dbc = @(x) (x+pi).*(x-2*pi); 
N.lbc = 0;

u = N \ 0; 

pass(j) = ( norm( u(:,d(3)) - N.dbc ) < tol); j = j + 1; 
pass(j) = ( norm( u(:,d(4)) - N.ubc ) < tol); j = j + 1;
pass(j) = ( norm( u(d(1),:) - N.lbc ) < tol); j = j + 1;
pass(j) = ( norm( u(d(2),:) - N.rbc ) < tol); j = j + 1;


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2), d); 
bdy = @(x,y) real(exp(x+1i*y)); 
N.lbc = @(y) bdy(d(1),y); N.rbc = @(y) bdy(d(2),y); 
N.dbc = @(x) bdy(x,d(3)); N.ubc = @(x) bdy(x,d(4)); 
u = N \ 0; exact = chebfun2(bdy, d);
pass(j) = ( norm( exact - u ) < tol ); j = j + 1; 

end
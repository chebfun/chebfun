function pass = test_domain( prefs ) 
%Check that boundary conditions are being imposed correctly. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.techPrefs.eps; 

d = [-2 2 -2 2];
N = chebop2(@(u) diff(u, 2, 1) + diff(u, 2, 2), d);

N.rbc = @(x) (2-x).*(2+x); 
N.ubc = 0; 
N.dbc = 0; 
N.lbc = 0;

u = N \ 0; 

pass(1) = ( norm( u(:,d(3)) - N.dbc.' ) < tol); 
pass(2) = ( norm( u(:,d(4)) - N.ubc.' ) < tol); 
pass(3) = ( norm( u(d(1),:) - N.lbc ) < tol); 
pass(4) = ( norm( u(d(2),:) - N.rbc ) < tol); 


d = [-pi 2*pi -2 5];
N = chebop2(@(u) diff(u, 2, 1) + diff(u, 2, 2), d);

N.rbc = 0; 
N.ubc = 0; 
N.dbc = @(x) (x+pi).*(x-2*pi); 
N.lbc = 0;

u = N \ 0; 

pass(5) = ( norm( u(:,d(3)) - N.dbc.' ) < tol);  
pass(6) = ( norm( u(:,d(4)) - N.ubc.' ) < 2*tol); 
pass(7) = ( norm( u(d(1),:) - N.lbc ) < 10*tol); 
pass(8) = ( norm( u(d(2),:) - N.rbc ) < 10*tol);


% Harmonic solution to the Laplace equation
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2), d); 
bdy = @(x,y) real(exp(x+1i*y)); 
N.lbc = @(y) bdy(d(1),y); N.rbc = @(y) bdy(d(2),y); 
N.dbc = @(x) bdy(x,d(3)); N.ubc = @(x) bdy(x,d(4)); 
u = N \ 0; exact = chebfun2(bdy, d);
pass(9) = ( norm( exact - u ) < 100*tol );

end
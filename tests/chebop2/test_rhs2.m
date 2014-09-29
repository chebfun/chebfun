function pass = test_rhs2( pref ) 

tol = 1e-14; 

% constant coefficients with rhs: 
op = @(u) laplacian(u);
N = chebop2( op ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
exact = N \ 1; 
% Other possible syntax: 
op = @(u) laplacian(u) - 1;
N = chebop2( op ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
u = N \ 0; 
pass(1) = norm( u - exact ) < tol; 

%%%  variable coefficient with rhs: 
% Standard syntax: 
N = chebop2( @(x,y,u) laplacian(u) ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
x = chebfun2(@(x,y) x ); 
exact = N \ sin(x);
% Other possible syntax
op = @(x,y,u) laplacian(u) - sin(x);
N = chebop2( op ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
u = N \ 0; 
pass(2) = norm( u - exact ) < tol;

%%%  Mixed rhs (part in operator, part on rhs). 
% Standard syntax: 
N = chebop2( @(x,y,u) laplacian(u) + sin(x.*(y + .1)) - 1) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
x = chebfun2(@(x,y) x ); 
exact = N \ 0;
% Other possible syntax
op = @(x,y,u) laplacian(u) + sin(x.*(y + .1));
N = chebop2( op ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
u = N \ 1; 
pass(3) = norm( u - exact ) < tol; 

end
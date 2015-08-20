function pass = test_bc( pref ) 
% Check that N.bc works 

tol = 1e-14; 

% constant coefficients with rhs: 
op = @(u) laplacian(u);
N = chebop2( op ) ; 
N.lbc = 0; N.rbc = 0; N.ubc = 0; N.dbc = 0; 
exact = N \ 1; 
% Other possible syntax: 
op = @(u) laplacian(u);
N = chebop2( op ) ; 
N.bc = 0;
u = N \ 1; 
pass(1) = norm( u - exact ) < tol; 

end
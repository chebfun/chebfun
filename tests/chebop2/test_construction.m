function pass = chebop2_construction
% Check that the chebop2 constructor is working correctly. 

tol = chebfun2pref('eps'); j = 1; 

% form laplacian: 
N = chebop2(@(u) laplacian(u)); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(j) = ( norm(u - 1) < tol ); j = j + 1;  

% second way to form laplacian:
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(j) = ( norm(u - 1) < tol ); j = j + 1; 

% on domains other than default 
N = chebop2(@(u) laplacian(u),[0 1 0 1]); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
myone = chebfun2(@(x,y) 1+0*x,[0 1 0 1]);   % check u is on correct domain.
pass(j) = ( norm(u - myone) < tol ); j = j + 1;  

% different boundary syntax 
N = chebop2(@(u) laplacian(u)); 
N.lbc = @(x) 1+0*x; N.rbc = chebfun(@(x) 1+0*x); 
N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(j) = ( norm(u - 1) < tol ); j = j + 1;

% more boundary syntax, just check it works. 
N = chebop2(@(u) laplacian(u));
N.lbc = @sin; N.rbc = chebfun(@(x) 1+0*x); 
N.ubc = 1; N.dbc = 1; 

% Neumann conditions. 
N = chebop2(@(u) laplacian(u));
N.lbc = @(x,u) diff(u) ; N.rbc = chebfun(@(x) 1+0*x); 
N.ubc = 1; N.dbc = 1; 
u = N \ 0 ; 

% Two ways to construct a chebop2 object. 
N1 = chebop2(@(x,y,u) diff(u,2,1) + diff(u,2,2) + x.*u); 
x = chebfun2(@(x,y) x); 
N2 = chebop2(@(u) diff(u,2,1) + diff(u,2,2) + x.*u);

C1 = N1.coeffs; C2 = N2.coeffs; 

end
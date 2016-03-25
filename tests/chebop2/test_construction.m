function pass = test_construction( prefs ) 
% Check that the chebop2 constructor is working correctly. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 10*prefs.techPrefs.chebfuneps;

% form laplacian: 
N = chebop2(@(u) laplacian(u)); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(1) = ( norm(u - 1) < tol );   

% second way to form laplacian:
N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(2) = ( norm(u - 1) < tol );  

% on domains other than default 
N = chebop2(@(u) laplacian(u),[0 1 0 1]); 
N.lbc = 1; N.rbc = 1; N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
myone = chebfun2(@(x,y) 1+0*x,[0 1 0 1]);   % check u is on correct domain.
pass(3) = ( norm(u - myone) < tol ); 

% different boundary syntax 
N = chebop2(@(u) laplacian(u)); 
N.lbc = @(x) 1+0*x; N.rbc = chebfun(@(x) 1+0*x); 
N.ubc = 1; N.dbc = 1; 
u = N \ 0; 
pass(4) = ( norm(u - 1) < tol ); 

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
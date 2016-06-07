function pass = test_neumann( pref )
% Check that we can impose Neumann boundary conditions. 
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 100*pref.cheb2Prefs.chebfun2eps;

d = [-2 3 -1 1];
N = chebop2(@(u) laplacian(u), d); 
N.lbc = @(y) y; 
N.rbc = @(y) 2*y; 
N.ubc = @(x) (x-d(1))/(d(2)-d(1))+1; 
N.dbc = @(x,u) diff(u)-((x-d(1))/(d(2)-d(1))+1); 

u = N \ 0;

x = chebfun('x', d(1:2)).';
uy = diff(u,1,1);
pass(1) = ( norm(u(:,d(4)) - N.ubc.') < tol ); 
pass(2) = ( norm(u(d(2),:) - N.rbc) < tol ); 
pass(3) = ( norm(uy(:,d(3))-((x-d(1))/(d(2)-d(1))+1)) < 5*tol );
pass(4) = ( norm(u(d(1),:) - N.lbc) < tol );
pass(5) = ( norm(laplacian(u)) < 400*tol );

% Nick Hale's example for Neumann conditions: 
N = chebop2(@(u) laplacian(u), [0 1 0 1]);
N.ubc = 0; N.rbc = 0; N.lbc = 0;
N.dbc = @(x, u) diff(u) - sin(2*pi*x);
u = N\0;
dudx = diffy(u);
bc = chebfun( @(x) sin(2*pi*x), [0,1] ); 
pass(6) = ( norm( dudx(:,0) - bc' ) < 10*tol ); 


%%
% d = [-2 3 -1 1];
% N = chebop2(@(u) lap(u), d); 
% N.dbc = @(y) y; 
% N.ubc = @(y) 2*y; 
% N.rbc = @(x) 3*(x-d(3))/(d(4)-d(3)) + 3; 
% N.lbc = @(x,u) diff(u)-(1*(x-d(3))/(d(4)-d(3))+1); 
% 
% u = N \ 0;
% 
% x = chebfun('x', d(3:4));
% ux = diff(u,1,2);
% pass(j) = ( norm(u(:,d(4)) - N.ubc) < tol ); j = j + 1;
% pass(j) = ( norm(u(:,d(3)) - N.dbc) < tol ); j = j + 1;
% pass(j) = ( norm(u(d(2),:) - N.rbc) < tol ); j = j + 1;
% pass(j) = ( norm(ux(d(1),:)-(1.5*(x-d(3))/(d(4)-d(3))+3)) < tol ); j = j + 1;
% pass(j) = ( norm(u(d(1),:) - N.lbc) < tol ); j = j + 1;
% pass(j) = ( norm(lap(u)) < 10*tol ); j = j + 1;


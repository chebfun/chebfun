function pass = chebop2_ADtest
% Test that the AD machinery is work correctly. 
% Alex Townsend, August 2013. 


tol = chebfun2pref('eps'); 
j = 1; 

% simple case. 
N = chebop2(@(u) diffx(u,2) + diffy(u,2));
pass(j) = ( norm( N.coeffs - [0 0 1; 0 0 0; 1 0 0]) < tol );
 j = j + 1;

% Helmoltz 
N = chebop2(@(u) diffx(u,2) + diffy(u,2) + pi*u);
pass(j) = ( norm( N.coeffs - [pi 0 1; 0 0 0; 1 0 0]) < tol );
 j = j + 1;

% Higher order 
N = chebop2(@(u) diffx(u,3) + diffx(diffy(u,1),2) + diffy(u,2) + pi*u);
pass(j) = ( norm( N.coeffs - [pi 0 1; 0 0 0; 0 1 0; 1 0 0]) < tol );
j = j + 1;

% Higher order 
N = chebop2(@(u) diffx(u,3) + diffx(diffy(u,1),2) + diffy(u,2) + pi*u);
pass(j) = ( norm( N.coeffs - [pi 0 1; 0 0 0; 0 1 0; 1 0 0]) < tol );
j = j + 1;

% Variable coefficients on zero derivatives.
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,1,2) + diff(u,2,1)  + x.*u);
CC = N.coeffs;
pass(j) = ( norm( CC{1,1} - x) < tol ); j = j + 1;
pass(j) = ( norm( CC{1,2} - 0*x - 1 ) < tol ); j = j + 1;
pass(j) = ( norm( CC{3,1} - 0*x - 1 ) < tol ); j = j + 1;


% % Variable coefficients on derivatives
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y);
N = chebop2(@(x,y,u) diff(u,1,2) + x.*diff(u,2,1));
CC = N.coeffs;
pass(j) = ( norm( CC{3,1} - x) < tol ); j = j + 1;
pass(j) = ( norm( CC{1,2} - 0*x - 1 ) < tol ); j = j + 1;
pass(j) = ( norm( CC{1,1} - 0*x ) < tol ); j = j + 1;




end
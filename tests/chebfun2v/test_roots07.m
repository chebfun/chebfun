function pass = test_roots07( pref ) 
% Check that the marching squares and Bezoutian agree with each other. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

%%
f = chebfun2(@(x,y)(x.^2+y.^2-1).*(x-1.1)); 
g = chebfun2(@(x,y)(25*x.*y-12).*(x-1.1)); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%% (Marching squares misses some solutions)
f = chebfun2(@(x,y)y.^4 + (-1)*y.^3 + (2*x.^2).*(y.^2) + (3*x.^2).*y + (x.^4)); 
g = @(x,y)y.^10-2*(x.^8).*(y.^2)+4*(x.^4).*y-2; 
g = chebfun2(@(x,y) g(2*x,2*(y+.5)));
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(f(r2(:,1),r2(:,2))) < tol ); j = j + 1; 
pass(j) = ( norm(g(r2(:,1),r2(:,2))) < 1000*tol ); j = j + 1;

%% (Marching squares misses a solution)
a=1e-9; rect = a*[-1 1 -1 1];
f = chebfun2(@(x,y)cos(x.*y/a^2)+sin(3*x.*y/a^2),rect); 
g = chebfun2(@(x,y)cos(y/a)-cos(2*x.*y/a^2),rect);
% r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
% pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
% pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;
pass(j) = ( norm(f(r2(:,1),r2(:,2))) < 1000*tol ); j = j + 1; 
pass(j) = ( norm(g(r2(:,1),r2(:,2))) < 1000*tol ); j = j + 1;

%% (Marching squares misses the roots of the edge)
% f = chebfun2(@(x,y)sin(3*pi*x).*cos(x.*y)); 
% g = chebfun2(@(x,y)sin(3*pi*y).*cos(sin(x.*y))); 
% r2 = roots([f;g]); 
% pass(j) = ( norm(f(r2(:,1),r2(:,2))) < tol ); j = j + 1; 
% pass(j) = ( norm(g(r2(:,1),r2(:,2))) < tol ); j = j + 1;


end
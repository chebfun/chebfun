function pass = test_roots09( pref ) 
% Check that the marching squares and Bezoutian agree with each other. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

%%
f = chebfun2(@(x,y)x.^2+y.^2-.9^2); 
g = chebfun2(@(x,y)sin(x.*y));
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%%
f = chebfun2(@(x,y)x.^2+y.^2-.49^2); 
g = chebfun2(@(x,y)(x-.1).*(x.*y-.2));
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

end

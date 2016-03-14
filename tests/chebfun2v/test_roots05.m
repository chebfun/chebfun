function pass = test_roots05( pref )
% Check that the marching squares and Bezoutian agree with each other. 
%%

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

%%
rect = 2*[-1 1 -1 1];
f = chebfun2(@(x,y)2*x.*y.*cos(y.^2).*cos(2*x)-cos(x.*y),rect); 
g = chebfun2(@(x,y)2*sin(x.*y.^2).*sin(3*x.*y)-sin(x.*y),rect); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

end
function pass = test_roots_slow( pref )
% Check that the marching squares and Bezoutian agree with each other. 
%%

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

%% (slow one)
f = chebfun2(@(x,y)exp(x-2*x.^2-y.^2).*sin(10*(x+y+x.*y.^2))); 
g = chebfun2(@(x,y)exp(-x+2*y.^2+x.*y.^2).*sin(10*(x-y-2*x.*y.^2))); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

end
function pass = test_roots03( pref )
% Check that the marching squares and Bezoutian agree with each other. 
%%

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

%% (Marching squares fails)
f = chebfun2(@(x,y)((x-.3).^2+2*(y+0.3).^2-1)); 
g = chebfun2(@(x,y)((x-.49).^2+(y+.5).^2-1).*((x+0.5).^2+(y+0.5).^2-1).*((x-1).^2+(y-0.5).^2-1)); 
%r1 = roots([f;g]); 
r2 = roots([f;g],'resultant'); 
pass(j) = ~( length(r2) - 4 ); j = j+1;
pass(j) = ( norm(f(r2(:,1),r2(:,2))) < tol ); j = j + 1; 
pass(j) = ( norm(g(r2(:,1),r2(:,2))) < 1e2*tol ); j = j + 1; 

%% (Marching Squares misses a root)
f = chebfun2(@(x,y)((x-0.1).^2+2*(y-0.1).^2-1).*((x+0.3).^2+2*(y-0.2).^2-1).*((x-0.3).^2+2*(y+0.15).^2-1).*((x-0.13).^2+2*(y+0.15).^2-1)); 
g = chebfun2(@(x,y)(2*(x+0.1).^2+(y+0.1).^2-1).*(2*(x+0.1).^2+(y-0.1).^2-1).*(2*(x-0.3).^2+(y-0.15).^2-1).*((x-0.21).^2+2*(y-0.15).^2-1)); 
% r1 = roots([f;g]); 
r2 = roots([f;g],'resultant'); 
pass(j) = ~( length(r2) - 45 ); j = j+1;
pass(j) = ( norm(f(r2(:,1),r2(:,2))) < 10*tol ); j = j + 1; 
pass(j) = ( norm(g(r2(:,1),r2(:,2))) < 100*tol ); j = j + 1; 

end
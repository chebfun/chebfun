function pass = test_roots01( pref )
% Check that the marching squares and Bezoutian agree with each other. 
% Uncomment tests if harder tests should be executed.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

%% 
f = chebfun2(@(x,y) 144*(x.^4+y.^4)-225*(x.^2+y.^2) + 350*x.^2.*y.^2+81); 
g = chebfun2(@(x,y) y-x.^6); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1; 

%% 
f = chebfun2(@(x,y)(y.^2-x.^3).*((y-0.7).^2-(x-0.3).^3).*((y+0.2).^2-(x+0.8).^3).*((y+0.2).^2-(x-0.8).^3)); 
g = chebfun2(@(x,y)((y+.4).^3-(x-.4).^2).*((y+.3).^3-(x-.3).^2).*((y-.5).^3-(x+.6).^2).*((y+0.3).^3-(2*x-0.8).^3)); 
r2 = roots([f;g],'resultant'); 
pass(j) = ~( length(r2) - 13 );  j = j + 1;
pass(j) = ( norm(f(r2(:,1),r2(:,2))) < 1e3*tol ); j = j + 1; 
pass(j) = ( norm(g(r2(:,1),r2(:,2))) < 1e3*tol ); j = j + 1;


%%
f = chebfun2(@(x,y)y.^2-x.^3); 
g = chebfun2(@(x,y)(y+.1).^3-(x-.1).^2); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%%
p = chebfun2(@(x,y) x - y + .5);
q = chebfun2(@(x,y) x + y );
r = roots([p; q]); 
pass(j) = norm( r - [-.25 .25] ) < tol; j = j + 1; 

p = chebfun2(@(x, y) y + x/2 + 1/10);
q = chebfun2(@(x, y) y - 2.1*x + 2);
r = roots([p ;  q], 'resultant');
pass(j) = norm(r - [0.730769230769231, -0.465384615384615]) < tol; j = j + 1;

end

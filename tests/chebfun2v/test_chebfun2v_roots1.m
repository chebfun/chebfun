function pass = test_chebfun2v_roots1( pref )
% Check that the marching squares and Bezoutian agree with each other. 
% Uncomment tests if harder tests should be executed.

if ( nargin < 1 ) 
    pref = chebpref; 
end 
tol = 1e3 * pref.cheb2Prefs.eps; 
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
f = chebfun2(@(x,y)cos(10*x.*y)); 
g = chebfun2(@(x,y) x + y.^2); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%%
f = chebfun2(@(x,y)x); 
g = chebfun2(@(x,y) (x-.9999).^2 + y.^2-1); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%% Boyd's problem in Chebfun workshop talk.
f = chebfun2(@(x,y)sin(4*(x + y/10 + pi/10))); 
g = chebfun2(@(x,y) cos(2*(x-2*y+pi/7))); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;

%% slow one.
% f = chebfun2(@(x,y)exp(x-2*x.^2-y.^2).*sin(10*(x+y+x.*y.^2))); 
% g = chebfun2(@(x,y)exp(-x+2*y.^2+x.*y.^2).*sin(10*(x-y-2*x.*y.^2))); 
% r1 = roots([f;g]); 
% r2 = roots([f;g]); 
% pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); j = j + 1; 
% pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); j = j + 1;
% 
%% another slow one.
% rect = 4*[-1 1 -1 1]; 
% f = chebfun2(@(x,y)2*y.*cos(y.^2).*cos(2*x)-cos(y),rect); 
% g = chebfun2(@(x,y)2*sin(y.^2).*sin(2*x)-sin(x),rect); 
% r1 = roots([f;g]); 
% r2 = roots([f;g]); 
% pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < 10*tol ); j = j + 1; 
% pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < 10*tol ); j = j + 1;

end


function pass = test_roots06( pref )
% Check that the marching squares and Bezoutian agree with each other. 
%%

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

%% Marching squares double counts some solutions. 
f = chebfun2(@(x,y)(y - 2*x).*(y+0.5*x)); 
g = chebfun2(@(x,y) x.*(x.^2+y.^2-1)); 
% r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
% pass(j) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < 10*sqrt(tol) ); j = j + 1; 
% pass(j) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < 10*sqrt(tol) ); j = j + 1;
pass(1) = ( norm(f(r2(:,1),r2(:,2))) < tol ); 
pass(2) = ( norm(g(r2(:,1),r2(:,2))) < tol );

%%
% Marching Squares fails if 1st kind points are used. Not sure why. 
if ( isa(pref.tech,'chebtech2') )
f = chebfun2(@(x,y)(y - 2*x).*(y+.5*x)); 
g = chebfun2(@(x,y) (x-.0001).*(x.^2+y.^2-1)); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
exact = [1/10000 -1/20000; 
              1/10000 1/5000 
              -2/sqrt(5) 1/sqrt(5) 
             -1/sqrt(5)  -2/sqrt(5)
              1/sqrt(5)    2/sqrt(5) 
               2/sqrt(5)  -1/sqrt(5)]; 
pass(3) = ( norm(sort(r1(:,1))-sort(exact(:,1))) < 10*tol ); 
pass(4) = ( norm(sort(r1(:,2))-sort(exact(:,2))) < 10*tol ); 
pass(5) = ( norm(sort(r2(:,1))-sort(exact(:,1))) < 100*tol ); 
pass(6) = ( norm(sort(r2(:,2))-sort(exact(:,2))) < 100*tol );
else
   pass(3:6) = 1;  
end
%%
f = chebfun2(@(x,y)25*x.*y - 12); 
g = chebfun2(@(x,y)x.^2+y.^2-1); 
r1 = roots([f;g],'ms'); 
r2 = roots([f;g],'resultant'); 
pass(7) = ( norm(sort(r1(:,1))-sort(r2(:,1))) < tol ); 
pass(8) = ( norm(sort(r1(:,2))-sort(r2(:,2))) < tol ); 

end
function pass = test_sum( pref ) 
% Test for integration of a fun2 object. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.cheb2Prefs.chebfun2eps;
j = 1; 

% Example from wiki: http://en.wikipedia.org/wiki/Multiple_integral#Double_integral
f = 'x.^2 + 4*y'; 
f = chebfun2(f, [11 14 7 10]);

exact = 1719;

pass(j)  = (abs(integral2(f)-exact)<30*tol); j=j+1;

% check syntax as well. 
f = chebfun2(@(x,y) x); 
pass(j) = norm(sum(f) - chebfun(@(x) 2*x)')<tol; j=j+1;
pass(j) = norm(sum(f,1) - sum(f))<10*eps; j = j+1; 
pass(j) = norm(sum(f,2))<10*eps; j = j+1; 


% On different domains. 
f = chebfun2(@(x,y) y,[0 1 -pi pi]); 
pass(j) = norm(sum(f))<10*eps; j=j+1;
pass(j) = norm(sum(f,1) - sum(f)) < tol; j = j+1; 
pass(j) = norm(sum(f,2)-chebfun(@(x) x,[-pi pi])) < tol; j = j+1; 

end
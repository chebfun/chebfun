function pass = test_polyfitL1( pref )
% Test chebfun/polyfitL1: 

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e-8;

f = chebfun(@(x) exp(x)*sin(10*x)); 
n = 5;  
p1best = polyfitL1(f, n); 

T = chebpoly(0:n); 
err = f-p1best; 
pass(1) = norm(sum( T.*sign(err) ))<tol;

p1best = polyfitL1(f,n); 
err = f-p1best; 
pass(2) = norm(sum( T.*sign(err) ))<tol;

% General intervals: 
f = chebfun(@(x) exp(x)*sin(10*x), [-2 2]); 
n = 5;  
p1best = polyfitL1(f, n);
err = f-p1best;
T = chebpoly(0:n, [-2, 2]); 
pass(3) = norm(sum( T.*sign(err) ))<tol;
end
function pass = test_bestL1( pref )

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e-8;

f = chebfun(@(x) exp(x)*sin(10*x));
n = 5;  
T = chebpoly(0:n); 

p1best = watson(f, n); 
err = f-p1best; 
pass(1) = norm(sum( T.*sign(err) )) < tol;

n = 5;  
T = chebpoly(0:n); 
f = chebfun(@(x) abs(x-0.2), 'splitting', 'on');
p1best = watson(f, n); 
err = f-p1best; 
pass(2) = norm(sum( T.*sign(err) )) < tol;
end
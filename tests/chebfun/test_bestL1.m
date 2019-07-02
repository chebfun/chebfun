function pass = test_bestL1( pref )

if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e-8;

<<<<<<< HEAD
f = chebfun(@(x) exp(x)*sin(10*x)); func = 'expsin10';
%f = chebfun(@(x) abs(x-0.2),'splitting','on'); func = 'abs02';
=======
f = chebfun(@(x) exp(x)*sin(10*x));
>>>>>>> f7e034a86a38a4293de15ec37dd59756c22da76d
n = 5;  
p1best = L1min_watson(f,n); 

T = chebpoly(0:n); 

<<<<<<< HEAD
err = f-p1best; 
% plot(err),shg
pass(1) = norm(sum( T.*sign(err) ))<1e-8

p1best = watson(f,n); 
err = f-p1best; 
pass(2) = norm(sum( T.*sign(err) ))<1e-8

end
=======
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
>>>>>>> f7e034a86a38a4293de15ec37dd59756c22da76d

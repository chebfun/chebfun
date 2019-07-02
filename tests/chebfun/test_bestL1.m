

function pass = test_bestL1(pref)

f = chebfun(@(x) exp(x)*sin(10*x)); func = 'expsin10';
%f = chebfun(@(x) abs(x-0.2),'splitting','on'); func = 'abs02';
n = 5;  
p1best = L1min_watson(f,n); 

T = chebpoly(0:n); 

err = f-p1best; 
% plot(err),shg
pass(1) = norm(sum( T.*sign(err) ))<1e-8

p1best = watson(f,n); 
err = f-p1best; 
pass(2) = norm(sum( T.*sign(err) ))<1e-8

end

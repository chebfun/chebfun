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

n = 10; 
p1best = polyfitL1(f, n); 
err = f-p1best; 
pass(2) = norm(sum( T.*sign(err) ))<tol;

% General intervals: 
f = chebfun(@(x) exp(x)*sin(10*x), [-2 2]); 
n = 5;  
p1best = polyfitL1(f, n);
err = f-p1best;
T = chebpoly(0:n, [-2, 2]); 
pass(3) = norm(sum( T.*sign(err) ))<tol;

% Check optimality in the easy case: 
f = chebfun(@(x) abs(x), 'splitting', 'on'); 
n = 7; 
p = polyfitL1(f, n);
err = f - p; 
T = chebpoly(0:n); 
pass(4) = norm( sum( T.*sign(err) ) ) < tol;

% General interval, easy case: 
dom = [-2 2]; 
f = chebfun(@(x) abs(x), dom, 'splitting', 'on'); 
n = 7; 
p = polyfitL1(f, n);
err = f - p; 
T = chebpoly(0:n,dom); 
pass(5) = norm( sum( T.*sign(err) ) ) < tol;

% Nonsymmetric general intervals: 
f = chebfun(@(x) sin(x-.1).*cos(3*x), [-3 2]); 
n = 11;  
T = chebpoly(0:n, [-3 2]); 
p1best = polyfitL1(f, n);
err = f-p1best;
pass(6) = norm(sum( T.*sign(err) ))<tol;

% Long domain: 
dom = [0 100];
f = chebfun(@(x) sin(x-.1)+cos(3*x), dom); 
n = 11;  
T = chebpoly(0:n, dom); 
p1best = polyfitL1(f, n);
err = f-p1best;
pass(7) = norm(sum( T.*sign(err) ))<tol;

% Large degree n: 
dom = [0 100];
f = chebfun(@(x) sin(x-.1)+cos(3*x), dom); 
n = 101;  
T = chebpoly(0:n, dom); 
p1best = polyfitL1(f, n);
err = f-p1best;
pass(8) = norm(sum( T.*sign(err) ))<tol;

end
function pass = test_diag( pref ) 
% Test diag

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.eps; 
j = 1; 

% checking diag
f = chebfun2(@(x,y) cos(x.*y)); 
d1 = chebfun(@(x) cos(x.^2)); 
d2 = chebfun(@(x) cos(x.*(x+.1)), [-1 .9]); 
pass(j) = norm( diag(f) - d1 ) < tol; j = j + 1; 
pass(j) = norm( diag(f, .1) - d2 ) < tol; j = j + 1; 

% Test trace: 
pass(j) = norm( sum(diag(f)) - trace(f) ) < tol; j = j + 1; 
pass(j) = norm( sum(d1) - trace(f) ) < tol; j = j + 1; 

% checking diag
f = chebfun2(@(x,y) cos(x.*y), [-3 3 -1 4]); 
d1 = chebfun(@(x) cos(x.^2), [-1 3]); 
d2 = chebfun(@(x) cos(x.*(x+.1)), [-1.1 3]); 
pass(j) = norm( diag(f) - d1 ) < tol; j = j + 1; 
pass(j) = norm( diag(f, .1) - d2 ) < tol; j = j + 1; 

% Test trace: 
pass(j) = norm( sum(diag(f)) - trace(f) ) < tol; j = j + 1; 
pass(j) = norm( sum(d1) - trace(f) ) < tol; j = j + 1; 

end
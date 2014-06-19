function pass = test_vertcat( prefs ) 
% Test vertical concatenation of CHEBFUN2 objects. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 10*prefs.eps;

f = chebfun2(@(x,y) x); 
g = chebfun2(@(x,y) y);
F = chebfun2v(@(x,y) x, @(x,y) y); 
pass(1) = norm( vertcat(f, g) - F ) < tol; 

d = [-1 0 -2 1];
f = chebfun2(@(x,y) x,d); 
g = chebfun2(@(x,y) y,d);
F = chebfun2v(@(x,y) x, @(x,y) y,d); 
pass(2) = norm( vertcat(f, g) - F ) < tol; 
pass(3) = norm( vertcat(f, g) - [ f ; g ] ) < tol; 

pass(4) = norm( vertcat(f) - f ) < tol; 

end
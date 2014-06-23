function pass = test_plus( pref ) 
% This tests the basic arithmetic operations on chebfun2 objects.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.eps; 
j = 1;

D = [-1 1 -1 1; -2 2 -2 2; -1 pi 0 2*pi];

for r = 1 : size(D,1)
    f = chebfun2(@(x,y) cos(x.*y), D(r,:));
    g = chebfun2(@(x,y) x + y + x.*y, D(r,:));
    
    FplusG = chebfun2(@(x,y) cos(x.*y) + (x + y + x.*y), D(r,:) );
    
    tolr = norm(D(r,:),inf)*tol;
    
    pass(j) = ( norm( f+g - FplusG ) < tolr ); j = j + 1;
end

% Check if chebfun2/plus compresses the rank: 
f = chebfun2(@(x,y) x); 
g = f + f; 
pass(4) = ( length(g) == length(f) ); 

% Check if chebfun2/plus compresses the rank: 
f = chebfun2(@(x,y) cos(x.*y)); 
g = f + f; 
pass(5) = ( length(g) == length(f) ); 

% Check adding a function with a small vscale works: 
f = chebfun2(@(x,y) 1e-100*x); 
g = f + f; 
pass(5) = ( length(g) == length(f) ); 
pass(6) = ( abs( vscale( g ) - 2e-100 ) < tol ); 

% Check adding a function with a large vscale works: 
f = chebfun2(@(x,y) 1e100*x); 
g = f + f; 
pass(7) = ( length(g) == length(f) ); 
pass(8) = ( abs( vscale( g ) - 2e100 )/2e100 < 2e-2 ); 

end
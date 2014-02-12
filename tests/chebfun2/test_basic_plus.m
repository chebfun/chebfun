function pass = test_basic_plus( pref ) 
% This tests the basic arithmetic operations on chebfun2 objects.

if ( nargin < 1 ) 
    pref = chebpref; 
end 
tol = 1e5 * pref.cheb2Prefs.eps; 
j = 1;

D = [-1 1 -1 1; -2 2 -2 2; -1 pi 0 2*pi];

for r = 1 : size(D,1)
    f = chebfun2(@(x,y) cos(x.*y), D(r,:));
    g = chebfun2(@(x,y) x + y + x.*y, D(r,:));
    
    FplusG = chebfun2(@(x,y) cos(x.*y) + (x + y + x.*y), D(r,:) );
    
    tolr = norm(D(r,:),inf)*tol;
    
    pass(j) = ( norm( f+g - FplusG ) < tolr ); j = j + 1;
end


end
function pass = test_uplus( pref ) 
% This tests the basic arithmetic operations on chebfun2 objects.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.eps; 
j = 1;

D = [-1 1 -1 1; -2 2 -2 2; -1 pi 0 2*pi];

for r = 1 : size(D,1)
    f = chebfun2(@(x,y) cos(x.*y), D(r,:));
    
    uplusF = f;
    tolr = norm(D(r,:),inf)*tol;
    
    pass(j) = ( norm( f - uplusF ) < tolr ); j = j + 1;
end


end
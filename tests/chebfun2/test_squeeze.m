function pass = test_squeeze( pref ) 
% Test SQUEEZE

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 10*pref.eps; 

% This does not squeeze: 
f = chebfun2(@(x,y) cos(x.*y)); 
g = squeeze( f );  
pass(1) = ( norm( f - g ) < 2*tol ); 

f = squeeze( chebfun2(@(x,y) cos(x)) ); 
g = chebfun(@(x) cos(x)); 
pass(2) = ( norm( f - g.' ) < 100*tol ); 

f = squeeze( chebfun2(@(x,y) cos(y), [-1 1 -2 3]) ); 
g = chebfun(@(x) cos(x), [-2 3]);
pass(3) = ( norm( f - g ) < tol ); 

end
function pass = test_chebfun2_times( pref ) 
% Test contour

if ( nargin == 0) 
    pref = chebfun2pref; 
end

tol = 1000*pref.eps; 
j = 1; 

f = chebfun2(@(x,y) cos(x.*y)); 
h = chebfun2(@(x,y) 2*cos(x.*y)); 
k = chebfun2(@(x,y) cos(x.*y).^2); 

pass(j) = norm( f.*2 - h ) < tol; j = j + 1; 
pass(j) = norm( f*2 - h ) < tol; j = j + 1; 
pass(j) = norm( 2*f - h ) < tol; j = j + 1; 
pass(j) = norm( 2.*f - h ) < tol; j = j + 1; 
pass(j) = norm( f.^2 - k ) < tol; j = j + 1; 
pass(j) = norm( f.*f - k ) < tol; j = j + 1;


end
function pass = test_complex( pref ) 
% Test complex

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = chebfun2(@(x,y) cos(x.*y)); 
pass(j) = norm( f - complex( f ) ) < tol; j = j + 1; 
pass(j) = norm( f + 1i*f - complex(f, f) ) < tol; j = j + 1; 

end
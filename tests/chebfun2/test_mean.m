function pass = test_mean( pref )
% Check the commands mean

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = chebfun2(@(x,y) sin((x-.1).*y)); 
pass(j) = ( abs( mean2( f ) - sum2(f)/4 ) < tol ); j = j + 1; 

f = chebfun2(@(x,y) sin((x-.1).*y), [-2 3 -1 4]); 
pass(j) = ( abs( mean2( f ) - sum2(f)/25 ) < tol ); j = j + 1; 


pass(j) = ( norm( mean( f , 1) - sum(f,1)/5 ) < tol ); j = j + 1; 
pass(j) = ( norm( mean( f ,2 ) - sum(f,2)/5 ) < tol ); j = j + 1; 

end
function pass = test_chebfun2_cdr( pref ) 
% Test CDR 

if ( nargin == 0) 
    pref = chebfun2pref; 
end

tol = 1000*pref.eps; 
j = 1; 

f = chebfun2(@(x,y) cos(x.*y));
[C, D, R] = cdr( f ); 
err = norm( f - C * D * R' );
pass(j) = err < tol; j = j + 1; 

f = chebfun2(@(x,y) cos(x.*y), [-3 4 -1 10]);
[C, D, R] = cdr( f ); 
err = norm( f - C * D * R' );
pass(j) = err < tol; j = j + 1; 

d = cdr( f ); 
pass(j) = norm( diag( D ) - d(:) ) < tol; j = j + 1; 

end
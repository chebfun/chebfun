function pass = test_chebfun_lu( pref ) 
% Test Chebfun LU command. 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = pref.eps; 
j = 1;

% Check accuracy on [-1,1]
x = chebfun(@(x) x); 
A = [1 x x.^2 x.^3 x.^4 x.^5]; 
[L, U, p] = lu( A ); 

pass(j) = norm( triu(U) - U ) < tol; j = j + 1; 
pass(j) = norm( A - L * U ) < tol; j = j + 1; 
pass(j) = norm( diag( L(p,:) ) - ones(size(L,2),1) ) < 10*tol; j = j + 1; 
pass(j) = norm( tril(L(p,:)) - L(p,:)) < 10*tol; j = j + 1; 

% Check accuracy on [-2,3]
x = chebfun(@(x) x, [-2 3]); 
A = [1 x x.^2 x.^3 x.^4 x.^5]; 
[L, U, p] = lu( A ); 

pass(j) = norm( triu(U) - U ) < tol; j = j + 1; 
pass(j) = norm( A - L * U ) < 100*tol; j = j + 1; 
pass(j) = norm( diag( L(p,:) ) - ones(size(L,2),1) ) < 10*tol; j = j + 1; 
pass(j) = norm( tril(L(p,:)) - L(p,:)) < 20*tol; j = j + 1; 

end
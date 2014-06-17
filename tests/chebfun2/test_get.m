function pass = test_get( pref ) 
% Test GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.eps; 
j = 1; 

f = chebfun2(@(x,y) cos(x.*y)); 
[C, D, R] = cdr( f ); 
pass(j) = norm( [-1 1 -1 1] - f.domain ) < tol; j = j + 1; 
pass(j) = norm( C - f.cols ) < tol; j = j + 1; 
pass(j) = norm( R - f.rows ) < tol; j = j + 1; 
pass(j) = norm( 1./diag(D)' - f.pivotValues ) < tol; j = j + 1; 
pass(j) = all( size(f.pivotLocations) == [length(f), 2] ); j = j + 1;

end
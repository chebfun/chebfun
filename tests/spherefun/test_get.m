function pass = test_get( pref ) 
% Test GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = spherefun(@(x,y,z) 1 + sin(pi*x.*y) + sin(pi*x.*z));
[C, D, R] = cdr( f ); 
pass(j) = norm( [-pi pi 0 pi] - f.domain ) < tol; j = j + 1; 
pass(j) = norm( C - f.cols ) < tol; j = j + 1; 
pass(j) = norm( R - f.rows ) < tol; j = j + 1; 
pass(j) = norm( 1./diag(D) - f.pivotValues ) < tol; j = j + 1; 
pass(j) = all( size(f.pivotLocations) == [length(f), 2] ); j = j + 1;
pass(j) = ( f.nonZeroPoles ) & ( abs(f(0,0,1)) > tol );

end
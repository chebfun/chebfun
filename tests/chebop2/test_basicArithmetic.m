function pass = test_basicArithmetic( pref )
% Check that we can do the basic arithmetic on chebop2 objects. 
% Alex Townsend, August 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = 10*pref.cheb2Prefs.chebfun2eps; 

N1 = chebop2(@(x,y,u) diff(u,2,2)); 
N2 = chebop2(@(x,y,u) x.*diff(u,2,1));
N = N1 + N2; 

EXACT = chebop2(@(x,y,u) diff(u,2,2) + x.*diff(u,2,1));
ERROR = N - EXACT;
err = 0; C = ERROR.coeffs;
for jj = 1:size(C,1)
    for kk = 1:size(C,2)
        err = err + norm(C{jj,kk});
    end
end
pass(1) = ( abs(err) < tol ); 



N1 = chebop2(@(u) diff(u,2,2)); 
N2 = chebop2(@(x,y,u) x.*diff(u,2,1));
N = N1 + N2; 

EXACT = chebop2(@(x,y,u) diff(u,2,2) + x.*diff(u,2,1));
ERROR = N - EXACT;
err = 0; C = ERROR.coeffs;
for jj = 1:size(C,1)
    for kk = 1:size(C,2)
        err = err + norm(C{jj,kk});
    end
end
pass(2) = ( abs(err) < tol );



N1 = chebop2(@(x,y,u) y.*diff(u,2,2)); 
N2 = chebop2(@(u) diff(u,2,1));
N = N1 + N2; 

EXACT = chebop2(@(x,y,u) y.*diff(u,2,2) + diff(u,2,1));
ERROR = N - EXACT;
err = 0; C = ERROR.coeffs;
for jj = 1:size(C,1)
    for kk = 1:size(C,2)
        err = err + norm(C{jj,kk});
    end
end
pass(3) = ( abs(err) < tol );


end
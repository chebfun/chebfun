function pass = chebop2_basicArithmetic
% Check that we can do the basic arithmetic on chebop2 objects. 
% Alex Townsend, August 2013. 

j = 1; 
tol = chebfun2pref('eps'); 


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
pass(j) = ( abs(err) < tol ); j = j+1; 



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
pass(j) = ( abs(err) < tol ); j = j+1; 



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
pass(j) = ( abs(err) < tol ); j = j+1; 


end
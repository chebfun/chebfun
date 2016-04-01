function pass = test_coeffs2( ) 
% test spherefun coeff2() command

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(lambda,theta) cos(theta)); 
c = [.5 0 .5]'; 
pass(1) = norm( coeffs2(f) - c ) < tol; 

f = spherefun(@(lambda,theta) sin(theta).*sin(lambda)); 
c = 1/4*[-1 0 1 ; 0 0 0 ; 1 0 -1]; 
pass(2) = norm( coeffs2(f) - c ) < tol;

f = spherefun(@(lambda,theta) sin(theta).*sin(lambda)+cos(theta)); 
c = 1/4*[-1 2 1 ; 0 0 0 ; 1 2 -1]; 
pass(3) = norm( coeffs2(f) - c ) < tol; 

end
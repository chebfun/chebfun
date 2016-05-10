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

f = spherefun(@(lambda,theta) sin(theta).^2.*sin(lambda).*cos(lambda)+cos(theta).^2); 
c = 1/16*[-1i 0 4 0 1i; 0 0 0 0 0; 2i 0 8 0 -2i; 0 0 0 0 0; -1i 0 4 0 1i]; 
pass(4) = norm( coeffs2(f) - c ) < tol;

% Construction from zeros matrix should maintain a zeros coefficient
% matrix of the same size.
f = spherefun(zeros(5,4));
pass(5) = norm( coeffs2(f) - zeros(5,4) ) == 0;


end
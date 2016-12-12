function pass = test_harmonic( ) 
% simple tests for cylindrical harmonic command.

tol = 3e2*chebfunpref().cheb2Prefs.chebfun2eps;


%example: build via jbessel w/known bessel root 
 %L = 5, m = 3; 
 F = diskfun(@(t,r) sqrt(2)/ (sqrt(pi)*abs(besselj(5+1, 15.7001740797116)))...
     *besselj(5, 15.7001740797116*r).*cos(5*t), 'polar'); 
 G = diskfun.harmonic(5,3); 
 pass(1) = norm( F-G ) < tol;
 
 F = diskfun(@(t,r) sqrt(2)/ (sqrt(pi)*abs(besselj(5+1, 15.7001740797116)))...
     *besselj(5, 15.7001740797116*r).*sin(5*t), 'polar'); 
 G = diskfun.harmonic(-5,3); 
 pass(2) = norm( F-G ) < tol;
 
%L = 0, m = 5; 
 F = diskfun(@(t,r) sqrt(2)/ (sqrt(2*pi)*abs(besselj(0+1, 14.9309177084877)))...
     *besselj(0, 14.9309177084877*r), 'polar'); 
 G = diskfun.harmonic(0,5); 
 pass(3) = norm(F-G) < tol;

%test orthonormality 

A=diskfun.harmonic(29,10);
B=diskfun.harmonic(6,23);
C = diskfun.harmonic(-4, 3); 
pass(4) = abs(sum2(A.*A)-1) < tol; 
pass(5) = abs(sum2(B.*B)-1) < tol; 
pass(6) = abs(sum2(A.*B)) < tol; 
pass(7) = abs(sum2(A.*C)) < tol; 
pass(8) = abs(sum2(C.*C)-1) < tol; 

%neumann conditions, test orthonormal
A = diskfun.harmonic(10,8, 'neumann'); 
B = diskfun.harmonic(-4, 3, 'neumann'); 
pass(8) = abs(sum2(A.*A)-1) < tol; 
pass(9) = abs(sum2(B.*B)-1) < tol; 
pass(10) = abs(sum2(A.*B)) < tol; 


end
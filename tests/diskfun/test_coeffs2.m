function pass = test_coeffs2( ) 
% test diskfun coeff2() command

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(t,r) r.^2, 'polar');  
c = [.5 0 .5]'; 
pass(1) = norm( coeffs2(f) - c ) < tol; 

f = diskfun(@(t,r) r.*sin(t), 'polar'); 
c = 1i/(2)*[0 0 0 ; 1 0 -1 ]; 
pass(2) = norm( coeffs2(f) - c ) < tol;


f = diskfun(@(t,r) r.^2+ r.*(4*cos(t)+2*sin(t)), 'polar'); 
c = 1/2*[ 0 1 0; 4+2*1i 0 4-2*1i; 0 1 0];
pass(3) = norm( coeffs2(f) - c ) < tol; 

f = diskfun(@(t,r) r.^3.*cos(3*t)+r.^2.*(sin(t)).^2, 'polar'); 
c = (1/8)*[0 -1 0 2 0 -1 0; 3 0 0 0 0 0 3; 0 -1 0 2 0 -1 0;1 0 0 0 0 0 1];
pass(4) = norm( coeffs2(f) - c ) < tol;

% Construction from zeros matrix should maintain a zeros coefficient
% matrix of the same size.
f = diskfun(zeros(5,4));
pass(5) = norm( coeffs2(f) - zeros(5,4) ) == 0;


end
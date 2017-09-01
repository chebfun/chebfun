function pass = test_coeffs2diskfun( ) 
% test diskfun coeff2diskfun() command

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun.coeffs2diskfun(0);
pass(1) = iszero(f);

f = diskfun(@(t,r) r.^2, 'polar');
g = diskfun.coeffs2diskfun(coeffs2(f));
pass(2) = norm(f-g) < tol;

f = diskfun(@(t,r) r.*sin(t), 'polar'); 
g = diskfun.coeffs2diskfun(coeffs2(f));
pass(3) = norm(f-g) < tol;

f = diskfun(@(x,y) exp(-10*((x-0.5/sqrt(2)).^2 + (y-0.5/sqrt(2)).^2)));
g = diskfun.coeffs2diskfun(coeffs2(f));
pass(4) = norm(f-g);

f = diskfun(@(t,r) r.*sin(t), 'polar'); 
c = 1i/(2)*[0 0 0 ; 1 0 -1 ]; 
pass(5) = norm( f - diskfun.coeffs2diskfun(c) ) < tol;

f = diskfun(@(t,r) r.^3.*cos(3*t)+r.^2.*(sin(t)).^2, 'polar'); 
c = (1/8)*[0 -1 0 2 0 -1 0; 3 0 0 0 0 0 3; 0 -1 0 2 0 -1 0;1 0 0 0 0 0 1];
pass(6) = norm( f - diskfun.coeffs2diskfun(c) ) < tol;

end
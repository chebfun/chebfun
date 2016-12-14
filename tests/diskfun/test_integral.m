function pass = test_integral( ) 
% Test INTEGRAL()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

%test empty diskfun
f = diskfun();
pass(1) = (integral(f) == 0); 
pass(2) = (integral(f, chebfun(@(x) x)) == 0);

%test unitcircle feature
f = diskfun(@(x,y) sin(x.^2+3*y));
pass(3) = ( abs(sum(f(:, 1))-integral(f, 'unitcircle')) < tol);
pass(4) = ( abs(integral(diskfun.harmonic(3,2, 'neumann'), 'unitcircle')) < tol);

% test integral along a contour (need more rigorous tests when complex-valued
% diskfuns are supported)
z = chebfun(@(x) .5*exp(1i*pi*x));
g = diskfun(@(x, y) exp(-2*x.^2-2*y.^2)); 
pass(5) = ( abs(.5*sum(g(:,.5))-integral(g, z)) < tol); 
z = chebfun(@(x) x*exp(1i*pi/4));
pass(6) = ( abs(sum(g(pi/4, :))-integral(g, z) )  < tol );
end
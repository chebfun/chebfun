function pass = test_integral(pref)
% Test CHEBFUN2V/INTEGRAL().

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

F = chebfun2v(@(x,y) -y, @(x,y) x);
c = chebfun(@(t) exp(1i*t), [-pi, pi]);
I = integral(F, c);
pass(1) = ( abs(I - 2*pi) < tol );

F = chebfun2v(@(x,y) x.*y, @(x,y) y);
c = chebfun(@(t) exp(1i*t), [0, pi/2]);
I = integral(F, c);
pass(2) = ( abs(I - 1/6) < tol );

end

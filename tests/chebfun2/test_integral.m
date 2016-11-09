function pass = test_integral(pref)
% Test INTEGRAL()

if ( nargin == 0)
    pref = chebfunpref;
end
tol = 100*pref.cheb2Prefs.chebfun2eps;

% circle
c = chebfun(@(t) exp(1i*t), [0 2*pi]);

f = chebfun2(@(x,y) 1+0*x);
pass(1) = abs( integral(f, c) - 2*pi ) < tol;

f = chebfun2(@(x,y) cos(x));
pass(2) = abs( integral(f, c) - 2*pi*besselj(0,1) ) < tol;

f = chebfun2(@(x,y) cos(x.*y));
exact = sum(chebfun(@(t) cos( cos(t).*sin(t) ), [0 2*pi]));
pass(3) = abs( integral(f, c) - exact ) < tol;

% Test integral along a real curve:
f = chebfun2(@(x,y) cos(x) + sin(y));
c = chebfun(@(t) t);
exact = sum(chebfun(@(t) cos(t)));
pass(4) = ( abs(integral(f, c) - exact) < tol );

end
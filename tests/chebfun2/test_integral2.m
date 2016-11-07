function pass = test_integral2(pref)
% Test INTEGRAL2()

if ( nargin == 0)
    pref = chebfunpref;
end
tol = 100*pref.cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) x.*y);
pass(1) = ( abs(integral2(f) - 0) < tol );

% Integrate over a smaller domain:
pass(2) = ( abs(integral2(f, [0, 1, 0, 1]) - 0.25) < tol );

f = chebfun2(@(x,y) x.^2 .* cos(y), [0, 3, -1, 1]);
Iexact = 9 * 2 * sin(1);
pass(3) = ( abs(integral2(f) - Iexact) < tol );

end
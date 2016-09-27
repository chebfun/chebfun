function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x);
pass(1) = ~isPeriodicTech(f);

f = chebfun3(@(x,y,z) cos(pi*x), 'trig');
pass(2) = isPeriodicTech(f);

% Test for complex function values:
f = chebfun3(@(x,y,z) exp(1i*pi*x), 'trig');
pass(3) = isPeriodicTech(f);

% On different domain:
f = chebfun3(@(x,y,z) cos(x).*sin(z), [ -pi, pi, -pi, pi, -pi, pi ], 'trig');
pass(4) = isPeriodicTech(f);

end

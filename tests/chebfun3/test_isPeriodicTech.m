function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) x);
pass(j) = ~isPeriodicTech(f);
j = j + 1;

f = chebfun3(@(x,y,z) cos(pi*x), 'trig');
pass(j) = isPeriodicTech(f);
j = j + 1;

% Test for complex function values:
f = chebfun3(@(x,y,z) exp(1i*pi*x), 'trig');
pass(j) = isPeriodicTech(f);
j = j + 1;

% On different domain:
f = chebfun3(@(x,y,z) cos(x).*sin(z), [ -pi, pi, -pi, pi, -pi, pi ], 'trig');
pass(j) = isPeriodicTech(f);

end

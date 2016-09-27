function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) cos(pi*x));
pass(1) = ~isPeriodicTech(f);

f = chebfun2(@(x,y) cos(pi*x), 'trig');
pass(2) = isPeriodicTech(f);

% Test on other domain
f = chebfun2(@(x,y) cos(x).*sin(y) + cos(3*x) + cos(5*x), ...
    [ -pi, pi, -pi, pi ], 'trig');
pass(3) = isPeriodicTech(f);

end

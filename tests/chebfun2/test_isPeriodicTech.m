function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

f = chebfun2(@(x,y) cos(pi*x));
pass(j) = ~isPeriodicTech(f);
j = j + 1;

f = chebfun2(@(x,y) cos(pi*x), 'trig');
pass(j) = isPeriodicTech(f);
j = j + 1;

% Test on other domain
f = chebfun2(@(x,y) cos(x).*sin(y) + cos(3*x) + cos(5*x), ...
    [ -pi, pi, -pi, pi ], 'trig');
pass(j) = isPeriodicTech(f);

end

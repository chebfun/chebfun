function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = chebfun2v(@(x,y) x, @(x,y) y);
pass(1) = ~isPeriodicTech(f);

f1 = chebfun2(@(x,y) cos(pi*x), [ -1, 1, -pi, pi ], 'trig');
f2 = chebfun2(@(x,y) sin(y), [ -1, 1, -pi, pi ], 'trig');
f = [f1; f2];
pass(2) = isPeriodicTech(f);
f = [f1; f1; f2];
pass(3) = isPeriodicTech(f);

end

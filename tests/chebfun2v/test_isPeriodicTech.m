function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

f = chebfun2v(@(x,y) x, @(x,y) y);
pass(j) = ~isPeriodicTech(f);
j = j + 1;

f1 = chebfun2(@(x,y) cos(pi*x), [ -1, 1, -pi, pi ], 'trig');
f2 = chebfun2(@(x,y) sin(y), [ -1, 1, -pi, pi ], 'trig');
f = [f1; f2];
pass(j) = isPeriodicTech(f);
j = j + 1;
f = [f1; f1; f2];
pass(j) = isPeriodicTech(f);

end

function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3v(@(x,y,z) x, @(x,y,z) y);
pass(1) = ~isPeriodicTech(f);

f1 = chebfun3(@(x,y,z) cos(pi*x), 'trig');
f2 = chebfun3(@(x,y,z) sin(pi*y), 'trig');
pass(2) = isPeriodicTech([f1; f2]);
pass(3) = isPeriodicTech([f1; f2; f2]);
F = [f1; f2];
pass(4) = isPeriodicTech([F; f2]);

end
function pass = test_isPeriodicTech(pref)
% Test isPeriodicTech().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3v(@(x,y,z) x, @(x,y,z) y);
pass(j) = ~isPeriodicTech(f);
j = j + 1;

f1 = chebfun3(@(x,y,z) cos(pi*x), 'trig');
f2 = chebfun3(@(x,y,z) sin(pi*y), 'trig');
pass(j) = isPeriodicTech([f1; f2]);
j = j + 1;
pass(j) = isPeriodicTech([f1; f2; f2]);
j = j + 1;
F = [f1; f2];
pass(j) = isPeriodicTech([F; f2]);

end
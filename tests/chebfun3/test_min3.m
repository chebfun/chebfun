function pass = test_min3(pref)
% Test @chebfun3/min3 command.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

a = 0.3;
b = -0.4322;
c = -0.83343;
f = chebfun3(@(x,y,z) (x-a).^2 + (y-b).^2 + (z-c).^2);
exactVal = 0;
[val, ~] = min3(f);
pass = abs(val - exactVal) < tol;

end
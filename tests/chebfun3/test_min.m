function pass = test_min(pref)
% Test the chebfun3/min command. 

if ( nargin < 1 ) 
    pref = chebfunpref;
end
tol = 1e11*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2+y.^2+z.^2);
exactMin = chebfun2(@(y,z) y.^2+z.^2);

minF = min(f); 
minF1 = min(f, []);
minF2 = min(f, [], 2);
minF3 = min(f, [], 3);
h = min(f, [], 4); % Should be f itself.

pass(1) = norm(minF - exactMin) < tol;
pass(2) = norm(minF1 - exactMin) < tol;
pass(3) = norm(minF2 - exactMin) < tol;
pass(4) = norm(minF3 - exactMin) < tol;
pass(5) = norm(h - f) < tol;

end
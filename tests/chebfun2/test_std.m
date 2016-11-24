function pass = test_std(pref)
% Test separableApprox/std command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

f1 = chebfun2(@(x,y) y.^2 + x); % f1 has y-std equal to sqrt(4/45).
f2 = chebfun2(@(x,y) x.^2 + y); % f2 has x-std equal to sqrt(4/45).

exact = chebfun(@(x) sqrt(4/45));

s1 = std(f1);
ss1 = std(f1, [], 1);
s2 = std(f2, [], 2);

pass(1) = norm(s1 - exact) < tol;
pass(2) = norm(ss1 - exact) < tol;
pass(3) = norm(s2 - exact) < tol;

end
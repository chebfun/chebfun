function pass = test_vscale(pref)

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

f = diskfunv(@(x,y) cos(x), @(x,y) 2*y);
pass(1) = (vscale(f) <= 2);

end

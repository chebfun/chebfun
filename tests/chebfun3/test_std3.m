function pass = test_std3(pref)
% Test the chebfun3/max command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(10); % (x,y,z)-std of any fixed value is zero.
exact = 0;
h = std3(f);
pass(1) = norm(h - exact) < tol;

end
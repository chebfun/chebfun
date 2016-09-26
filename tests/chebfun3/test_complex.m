function pass = test_complex(pref) 
% Test complex

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) sin(x.*y.*z));
g = chebfun3(@(x,y,z) cos(x.*y.*z));
h = chebfun3(@(x,y,z) sin(x.*y.*z) + 1i*cos(x.*y.*z));
pass = norm(h - complex(f, g)) < tol;

end
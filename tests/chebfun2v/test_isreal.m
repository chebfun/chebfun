function pass = test_isreal(pref)
% Test chebfun2v/isreal.

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = chebfun2v(@(x,y) x, @(x,y) y);
pass(1) = isreal(f);

f = chebfun2v(@(x,y) x, @(x,y) 1i*y);
pass(2) = ~isreal(f);

f = real(chebfun2v(@(x,y) 1i*x + y, @(x,y) y-x));
pass(3) = isreal(f);

end
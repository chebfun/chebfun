function pass = test_isreal(pref)
% Test chebfun2v/isreal.

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

f = chebfun2v(@(x,y) x, @(x,y) y);
pass(j) = isreal(f);
j = j+1;

f = chebfun2v(@(x,y) x, @(x,y) 1i*y);
pass(j) = ~isreal(f);
j = j+1;

f = real(chebfun2v(@(x,y) 1i*x + y, @(x,y) y-x));
pass(j) = isreal(f);

end
function pass = test_cumsum3(pref)

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Check cumsum on cube domain
[x, y, z] = cheb.xyz;
f = x;
g = cumsum3(f);
gExact = (y+1).*(z+1).*(x.^2/2 - 1/2);
pass(1) = norm(g - gExact) < tol;

end
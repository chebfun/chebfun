function pass = test_poldec(pref)

if ( nargin < 1 )
    pref = chebfunpref; 
end 
tol = 1000*pref.cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) exp(x.*y));
[U,H] = poldec(f);

% check that U is partial isometry 
pass(1) = sum( abs( svd(U)-1 ) ) < tol;
% check that residual is small
pass(2) = norm(U*H - f) < tol;


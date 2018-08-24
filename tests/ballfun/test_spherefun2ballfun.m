function pass = test_spherefun2ballfun( pref ) 
% Test with function cos(cos(lam)*sin(th))

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

f = spherefun(@(lam,th)cos(cos(lam).*sin(th)));
[p,n] = size(coeffs2(f));
exact = ballfun(@(r,lam,th)cos(cos(lam).*sin(th)),[2,n,p]);
g = ballfun.spherefun2ballfun(f);
pass(1) = norm( g - exact ) < tol;

end

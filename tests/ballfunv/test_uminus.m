function pass = test_uminus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
F1 = ballfun(@(r,lam,th)r.^2.*cos(lam));
F2 = ballfun(@(r,lam,th)r.*sin(lam).^3.*sin(th));
F3 = ballfun(@(r,lam,th)cos(lam).*sin(th));
F = ballfunv(F1,F2,F3);
G = -F;
pass(1) = norm( F + G ) < tol;
end

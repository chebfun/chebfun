function pass = test_uplus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function rand = +rand
F1 = ballfun(@(r,lam,th)r.*cos(lam));
F2 = ballfun(@(r,lam,th)r.*sin(lam).*sin(th));
F3 = ballfun(@(r,lam,th)cos(lam).*cos(th));
F = ballfunv(F1,F2,F3);
G = +F;
pass(1) = norm( F - G ) < tol;
end

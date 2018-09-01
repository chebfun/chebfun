function pass = test_uplus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function rand = +rand
S = [19,12,27];
F1 = ballfun(@(r,lam,th)r.*cos(lam),S);
F2 = ballfun(@(r,lam,th)r.*sin(lam).*sin(th),S);
F3 = ballfun(@(r,lam,th)cos(lam).*cos(th),S);
F = ballfunv(F1,F2,F3);
G = +F;
pass(1) = norm( F - G ) < tol;
end

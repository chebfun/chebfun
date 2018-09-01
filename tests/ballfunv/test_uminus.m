function pass = test_uminus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function rand : rand + -rand = 0
S = [20,21,22];
F1 = ballfun(@(r,lam,th)r.^2.*cos(lam),S);
F2 = ballfun(@(r,lam,th)r.*sin(lam).^3.*sin(th),S);
F3 = ballfun(@(r,lam,th)cos(lam).*sin(th),S);
F = ballfunv(F1,F2,F3);
G = -F;
pass(1) = norm( F + G ) < tol;
end

function pass = test_minus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
F1 = ballfun(@(r,lam,th)r.*cos(lam));
F2 = ballfun(@(r,lam,th)r.*sin(lam).*sin(th));
F3 = ballfun(@(r,lam,th)cos(lam).*cos(th));
f = ballfunv(F1,F2,F3);
g = 2*f;
pass(1) = norm(g-f - f)<tol;

% Example 2
f1 = ballfun(@(r,lam,th)r.*cos(th));
f2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam));
f3 = ballfun(@(r,lam,th)r.*cos(th).^2);
f = ballfunv(f1,f2,f3);
g = f-5;
e1 = ballfun(@(r,lam,th)r.*cos(th)-5);
e2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam)-5);
e3 = ballfun(@(r,lam,th)r.*cos(th).^2-5);
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_minus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

S = [20,21,22];

% Example 1
F1 = ballfun(@(r,lam,th)r.*cos(lam),S);
F2 = ballfun(@(r,lam,th)r.*sin(lam).*sin(th),S);
F3 = ballfun(@(r,lam,th)cos(lam).*cos(th),S);
f = ballfunv(F1,F2,F3);
g = 2*f;
pass(1) = norm(g-f - f)<tol;

% Example 2
f1 = ballfun(@(r,lam,th)r.*cos(th),S);
f2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam),S);
f3 = ballfun(@(r,lam,th)r.*cos(th).^2,S);
f = ballfunv(f1,f2,f3);
g = f-5;
e1 = ballfun(@(r,lam,th)r.*cos(th)-5,S);
e2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam)-5,S);
e3 = ballfun(@(r,lam,th)r.*cos(th).^2-5,S);
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

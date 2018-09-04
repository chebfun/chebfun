function pass = test_plus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f1 = ballfun(@(r,lam,th)r.*cos(th));
f2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam));
f3 = ballfun(@(r,lam,th)r.*cos(th).^2);
f = ballfunv(f1,f2,f3);
g = f+f;
e1 = ballfun(@(r,lam,th)2*r.*cos(th));
e2 = ballfun(@(r,lam,th)2*r.*cos(th).*sin(lam));
e3 = ballfun(@(r,lam,th)2*r.*cos(th).^2);
exact = ballfunv(e1,e2,e3);
pass(1) = norm(g-exact)<tol;

% Example 2
f1 = ballfun(@(r,lam,th)r.*cos(th));
f2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam));
f3 = ballfun(@(r,lam,th)r.*cos(th).^2);
f = ballfunv(f1,f2,f3);
g = f+1;
e1 = ballfun(@(r,lam,th)r.*cos(th)+1);
e2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam)+1);
e3 = ballfun(@(r,lam,th)r.*cos(th).^2+1);
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

% Example 3
f1 = ballfun(@(r,lam,th)r.*cos(th));
f2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam));
f3 = ballfun(@(r,lam,th)r.*cos(th).^2);
f = ballfunv(f1,f2,f3);
g = 3+f;
e1 = ballfun(@(r,lam,th)r.*cos(th)+3);
e2 = ballfun(@(r,lam,th)r.*cos(th).*sin(lam)+3);
e3 = ballfun(@(r,lam,th)r.*cos(th).^2+3);
exact = ballfunv(e1,e2,e3);
pass(3) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

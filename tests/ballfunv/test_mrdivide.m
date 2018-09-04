function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

zero = cheb.galleryballfun('zero');

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2);
V = ballfunv(f,zero,zero);
W = V/2;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/2);
exact = ballfunv(g,zero,zero);
pass(1) = norm(W-exact)<tol;

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2);
V = ballfunv(f,f,f);
W = V/(-3);
g = ballfun(@(r,lam,th)-r.^2.*cos(lam).*sin(th).^2/3);
exact = ballfunv(g,g,g);
pass(2) = norm(W-exact)<tol;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2);
V = ballfunv(f,zero,f);
W = V/1i;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/1i);
exact = ballfunv(g,zero,g);
pass(3) = norm(W-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

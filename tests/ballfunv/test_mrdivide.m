function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

S = [20,21,22];
zero = cheb.galleryballfun('zero',S);

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,zero,zero);
W = V/2;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/2,S);
exact = ballfunv(g,zero,zero);
pass(1) = norm(W-exact)<tol;

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,f,f);
W = V/(-3);
g = ballfun(@(r,lam,th)-r.^2.*cos(lam).*sin(th).^2/3,S);
exact = ballfunv(g,g,g);
pass(2) = norm(W-exact)<tol;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,zero,f);
W = V/1i;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/1i,S);
exact = ballfunv(g,zero,g);
pass(3) = norm(W-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

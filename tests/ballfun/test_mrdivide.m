function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [20,21,22];

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/2;
exact = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/2,S);
pass(1) = norm( g - exact ) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/(-3);
exact = ballfun(@(r,lam,th)-r.^2.*cos(lam).*sin(th).^2/3,S);
pass(2) = norm( g - exact ) < tol;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/1i;
exact = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/1i,S);
pass(3) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

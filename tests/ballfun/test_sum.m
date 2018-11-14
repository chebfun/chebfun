function pass = test_sum( pref ) 
% Test with function cos(cos(lam)*sin(th))

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

%% Integrate over r

% Example 1
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th), 'polar');
g = sum(f, 1);
exact = spherefun(@(lam,th)cos(lam).*sin(th)/2,'vectorize');
pass(1) = norm(g-exact) < tol;

% Example 2
f = ballfun(@(r,lam,th)1, 'polar');
g = sum(f, 1);
exact = spherefun(@(lam,th)1,'vectorize');
pass(2) = norm(g-exact) < tol;

%% Integrate over lambda

% Example 3
f = ballfun(@(r,lam,th)(r.*sin(lam).*sin(th)).^2, 'polar');
g = sum(f, 2);
exact = diskfun(@(th,r)pi*r.^2.*sin(th).^2,'polar','vectorize');
pass(3) = norm(g-exact) < tol;

%% Integrate over theta

% Example 4
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th), 'polar');
g = sum(f, 3);
exact = diskfun(@(lam,r)2*r.*cos(lam),'polar','vectorize');
pass(4) = norm(g-exact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_rotate( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),'polar');
g = rotate(f,pi,0,0);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam-pi),'polar');
pass(1) = norm(g-exact) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),'polar');
g = rotate(f,pi/2,0,0);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam-pi/2),'polar');
pass(2) = norm(g-exact) < tol;

% Example 3
f = ballfun(@(r,lam,th)r.*cos(th),'polar');
g = rotate(f,0,pi,0);
exact = ballfun(@(r,lam,th)-r.*cos(th),'polar');
pass(3) = norm(g-exact) < tol;

% Example 4
f = ballfun(@(x,y,z)y);
g = rotate(f,-pi,0,pi);
pass(4) = norm(f-g) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
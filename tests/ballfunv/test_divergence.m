function pass = test_divergence( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1 : (x,y,z)
Vx = ballfun(@(r,lam,th)r.*sin(th).*cos(lam), 'spherical');
Vy = ballfun(@(r,lam,th)r.*sin(th).*sin(lam), 'spherical');
Vz = ballfun(@(r,lam,th)r.*cos(th), 'spherical');
V = ballfunv(Vx,Vy,Vz);
f = divergence(V);
exact = ballfun(@(r,lam,th)3, 'spherical');
pass(1) = norm(f-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
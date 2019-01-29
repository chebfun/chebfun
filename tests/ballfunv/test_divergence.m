function pass = test_divergence( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1: (x,y,z)
Vx = ballfun(@(r,lam,th)r.*sin(th).*cos(lam), 'spherical');
Vy = ballfun(@(r,lam,th)r.*sin(th).*sin(lam), 'spherical');
Vz = ballfun(@(r,lam,th)r.*cos(th), 'spherical');
V = ballfunv(Vx,Vy,Vz);
f = divergence(V);
exact = ballfun(@(r,lam,th)3, 'spherical');
pass(1) = norm(f-exact)  <tol;

% Exemple 2: check divergence theorem
v = ballfunv(@(x,y,z)sin(x), @(x,y,z)x.*z, @(x,y,z)cos(z));
lhs = sum3(div(v));
rhs = sum2(dot(v(1,:,:, 'spherical'), spherefunv.unormal));
pass(2) = abs(lhs-rhs) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
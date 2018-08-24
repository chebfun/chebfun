function pass = test_diff( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1 : dr 
f = ballfun(@(r,lam,th)r);
g = diff(f,1);
exact = ballfun(@(r,lam,th)1+0*r);
pass(1) = norm(f - exact) < tol;

% Example 2 : dlambda 
f = ballfun(@(r,lam,th)cos(lam));
g = diff(f,2);
exact = ballfun(@(r,lam,th)-sin(lam));
pass(2) = norm(f - exact) < tol;

% Example 3 : dtheta 
f = ballfun(@(r,lam,th)sin(th));
g = diff(f,3);
exact = ballfun(@(r,lam,th)cos(th));
pass(3) = norm(f - exact) < tol;

% Example 4 : dr^2 
f = ballfun(@(r,lam,th)r.^3.*sin(th));
g = diff(f,1,2);
exact = ballfun(@(r,lam,th)6*r.*sin(th));
pass(4) = norm(f - exact) < tol;

% Example 5 : dth^2
f = ballfun(@(r,lam,th)cos(th));
g = diff(f,3,2);
exact = ballfun(@(r,lam,th)-cos(th));
pass(5) = norm(f - exact) < tol;

% Example 6 : dlam^3
f = ballfun(@(r,lam,th)sin(lam));
g = diff(f,2,3);
exact = ballfun(@(r,lam,th)-cos(lam));
pass(6) = norm(f - exact) < tol;

% Example 7 : dx
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam));
g = diff(f,1,"cart");
exact = ballfun(@(r,lam,th)1);
pass(7) = norm(f - exact) < tol;

% Example 8 : dx
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*cos(lam).^2);
g = diff(f,1,"cart");
exact = ballfun(@(r,lam,th)2*r.*sin(th).*cos(lam));
pass(8) = norm(f - exact) < tol;

% Example 9 : dy
f = ballfun(@(r,lam,th)r.*sin(th).*sin(lam));
g = diff(f,2,"cart");
exact = ballfun(@(r,lam,th)1);
pass(9) = norm(f - exact) < tol;

% Example 10 : dy
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*sin(lam).^2);
g = diff(f,2,"cart");
exact = ballfun(@(r,lam,th)2*r.*sin(th).*sin(lam));
pass(10) = norm(f - exact) < tol;

% Example 11 : dy^2
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*sin(lam).^2);
g = diff(f,2,2,"cart");
exact = ballfun(@(r,lam,th)2);
pass(11) = norm(f - exact) < tol;

% Example 12 : dz
f = ballfun(@(r,lam,th)r.*cos(th));
g = diff(f,3,"cart");
exact = ballfun(@(r,lam,th)1);
pass(12) = norm(f - exact) < tol;

% Example 13 : dz^2
f = ballfun(@(r,lam,th)r.^2.*cos(th).^2);
g = diff(f,3,2,"cart");
exact = ballfun(@(r,lam,th)2);
pass(13) = norm(f - exact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

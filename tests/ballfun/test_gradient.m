function pass = test_gradient( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1: xyz
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th).*r.*cos(th),'spherical');
v = gradient(f);
exactx = ballfun(@(r,lam,th)r.*sin(lam).*sin(th).*r.*cos(th),'spherical');
exacty = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*cos(th),'spherical');
exactz = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th),'spherical');
exact = ballfunv(exactx,exacty,exactz);
pass(1) = norm( v - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
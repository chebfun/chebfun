function pass = test_curl( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1 : (0,x,xy)
Vx = ballfun(@(x,y,z)0);
Vy = ballfun(@(r,lam,th)r.*sin(th).*cos(lam), 'spherical');
Vz = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam), 'spherical');
V = ballfunv(Vx,Vy,Vz);
W = curl(V);
Exactx = ballfun(@(r,lam,th)r.*sin(th).*cos(lam), 'spherical');
Exacty = ballfun(@(r,lam,th)-r.*sin(th).*sin(lam), 'spherical');
Exactz = ballfun(@(r,lam,th)1, 'spherical');
Exact = ballfunv(Exactx,Exacty,Exactz);
pass(1) = norm(W-Exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
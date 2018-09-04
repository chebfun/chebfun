function pass = test_curl( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1 : (0,x,xy)
Vx = cheb.galleryballfun("zero");
Vy = ballfun(@(r,lam,th)r.*sin(th).*cos(lam));
Vz = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam));
V = ballfunv(Vx,Vy,Vz);
W = curl(V);
Exactx = ballfun(@(r,lam,th)r.*sin(th).*cos(lam));
Exacty = ballfun(@(r,lam,th)-r.*sin(th).*sin(lam));
Exactz = ballfun(@(r,lam,th)1);
Exact = ballfunv(Exactx,Exacty,Exactz);
pass(1) = norm(W-Exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
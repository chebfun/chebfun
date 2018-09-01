function pass = test_ballfunv2spherical( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

S = [31,32,33];
zero = cheb.galleryballfun('zero',S);

% Example 1 : (x,0,0)
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
v = ballfunv(f,zero,zero);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*sin(th).*cos(lam),S);
Vlam = ballfun(@(r,lam,th)-r.*sin(th).*cos(lam).*sin(lam),S);
Vth = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*cos(th).*cos(lam),S);
pass(1) = norm(Wr-Vr)<tol;
pass(2) = norm(Wlam-Vlam)<tol;
pass(3) = norm(Wth-Vth)<tol;

% Example 2 : (0,y,0)
f = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),S);
v = ballfunv(zero,f,zero);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*sin(th).*sin(lam),S);
Vlam = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*cos(lam),S);
Vth = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*cos(th).*sin(lam),S);
pass(4) = norm(Wr-Vr)<tol;
pass(5) = norm(Wlam-Vlam)<tol;
pass(6) = norm(Wth-Vth)<tol;

% Example 3 : (0,0,z)
f = ballfun(@(r,lam,th)r.*cos(th),S);
v = ballfunv(zero,zero,f);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*cos(th).*cos(th),S);
Vlam = zero;
Vth = ballfun(@(r,lam,th)-r.*cos(th).*sin(th),S);
pass(7) = norm(Wr-Vr)<tol;
pass(8) = norm(Wlam-Vlam)<tol;
pass(9) = norm(Wth-Vth)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
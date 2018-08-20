function pass = test_ballfunv2spherical()

S = [31,32,33];
zero = cheb.galleryballfun('zero',S);

% Example 1 : (x,0,0)
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
v = ballfunv(f,zero,zero);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*sin(th).*cos(lam),S);
Vlam = ballfun(@(r,lam,th)-r.*sin(th).*cos(lam).*sin(lam),S);
Vth = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*cos(th).*cos(lam),S);
pass(1) = isequal(Wr,Vr) && isequal(Wlam,Vlam) && isequal(Wth,Vth);

% Example 2 : (0,y,0)
f = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),S);
v = ballfunv(zero,f,zero);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*sin(th).*sin(lam),S);
Vlam = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*cos(lam),S);
Vth = ballfun(@(r,lam,th)r.*sin(th).*sin(lam).*cos(th).*sin(lam),S);
pass(2) = isequal(Wr,Vr) && isequal(Wlam,Vlam) && isequal(Wth,Vth);

% Example 3 : (0,0,z)
f = ballfun(@(r,lam,th)r.*cos(th),S);
v = ballfunv(zero,zero,f);
[Wr,Wlam,Wth] = ballfunv2spherical(v);
Vr = ballfun(@(r,lam,th)r.*cos(th).*cos(th),S);
Vlam = zero;
Vth = ballfun(@(r,lam,th)-r.*cos(th).*sin(th),S);
pass(3) = isequal(Wr,Vr) && isequal(Wlam,Vlam) && isequal(Wth,Vth);

if (nargout > 0)
    pass = all(pass(:));
end
end
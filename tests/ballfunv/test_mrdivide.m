function pass = test_mrdivide( ) 
S = [20,21,22];
zero = cheb.galleryballfun('zero',S);

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,zero,zero);
W = V/2;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/2,S);
exact = ballfunv(g,zero,zero);
pass(1) = isequal(W,exact);

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,f,f);
W = V/(-3);
g = ballfun(@(r,lam,th)-r.^2.*cos(lam).*sin(th).^2/3,S);
exact = ballfunv(g,g,g);
pass(2) = isequal(W,exact);

% Example 3
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
V = ballfunv(f,zero,f);
W = V/1i;
g = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/1i,S);
exact = ballfunv(g,zero,g);
pass(3) = isequal(W,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

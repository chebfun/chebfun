function pass = test_diff( )
S=[39,40,41];

% Example 1 : dr 
f = ballfun(@(r,lam,th)r,S);
g = diff(f,1);
exact = ballfun(@(r,lam,th)1+0*r,S);
pass(1) = isequal(g,exact);

% Example 2 : dlambda 
f = ballfun(@(r,lam,th)cos(lam),S);
g = diff(f,2);
exact = ballfun(@(r,lam,th)-sin(lam),S);
pass(2) = isequal(g,exact);

% Example 3 : dtheta 
f = ballfun(@(r,lam,th)sin(th),S);
g = diff(f,3);
exact = ballfun(@(r,lam,th)cos(th),S);
pass(3) = isequal(g,exact);

% Example 4 : dr^2 
f = ballfun(@(r,lam,th)r.^3.*sin(th),S);
g = diff(f,1,2);
exact = ballfun(@(r,lam,th)6*r.*sin(th),S);
pass(4) = isequal(g,exact);

% Example 5 : dth^2
f = ballfun(@(r,lam,th)cos(th),S);
g = diff(f,3,2);
exact = ballfun(@(r,lam,th)-cos(th),S);
pass(5) = isequal(g,exact);

% Example 6 : dlam^3
f = ballfun(@(r,lam,th)sin(lam),S);
g = diff(f,2,3);
exact = ballfun(@(r,lam,th)-cos(lam),S);
pass(6) = isequal(g,exact);

% Example 7 : dx
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
g = diff(f,1,"cart");
exact = ballfun(@(r,lam,th)1,S);
pass(7) = isequal(g,exact);

% Example 8 : dx
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*cos(lam).^2,S);
g = diff(f,1,"cart");
exact = ballfun(@(r,lam,th)2*r.*sin(th).*cos(lam),S);
pass(8) = isequal(g,exact);

% Example 9 : dy
f = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),S);
g = diff(f,2,"cart");
exact = ballfun(@(r,lam,th)1,S);
pass(9) = isequal(g,exact);

% Example 10 : dy
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*sin(lam).^2,S);
g = diff(f,2,"cart");
exact = ballfun(@(r,lam,th)2*r.*sin(th).*sin(lam),S);
pass(10) = isequal(g,exact);

% Example 11 : dy^2
f = ballfun(@(r,lam,th)r.^2.*sin(th).^2.*sin(lam).^2,S);
g = diff(f,2,2,"cart");
exact = ballfun(@(r,lam,th)2,S);
pass(11) = isequal(g,exact);

% Example 12 : dz
f = ballfun(@(r,lam,th)r.*cos(th),S);
g = diff(f,3,"cart");
exact = ballfun(@(r,lam,th)1,S);
pass(12) = isequal(g,exact);

% Example 13 : dz^2
f = ballfun(@(r,lam,th)r.^2.*cos(th).^2,S);
g = diff(f,3,2,"cart");
exact = ballfun(@(r,lam,th)2,S);
pass(13) = isequal(g,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

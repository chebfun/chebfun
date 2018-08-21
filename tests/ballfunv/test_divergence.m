function pass = test_divergence()

S = [23,41,37];

% Example 1 : (x,y,z)
Vx = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
Vy = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),S);
Vz = ballfun(@(r,lam,th)r.*cos(th),S);
V = ballfunv(Vx,Vy,Vz);
f = divergence(V);
exact = ballfun(@(r,lam,th)3,S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
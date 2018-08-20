function pass = test_grad()

% Example 1: xyz
S = [38,37,40];
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th).*r.*cos(th),S);
v = grad(f);
exactx = ballfun(@(r,lam,th)r.*sin(lam).*sin(th).*r.*cos(th),S);
exacty = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*cos(th),S);
exactz = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th),S);
exact = ballfunv(exactx,exacty,exactz);
pass(1) = isequal(v,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
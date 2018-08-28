function pass = test_rotate()

% Example 1
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam));
g = rotate(f,pi,0,0);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam+pi));
end
function pass = test_mrdivide( ) 
S = [20,21,22];

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/2;
exact = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/2,S);
pass(1) = isequal(g,exact);

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/(-3);
exact = ballfun(@(r,lam,th)-r.^2.*cos(lam).*sin(th).^2/3,S);
pass(2) = isequal(g,exact);

% Example 3
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = f/1i;
exact = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2/1i,S);
pass(3) = isequal(g,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

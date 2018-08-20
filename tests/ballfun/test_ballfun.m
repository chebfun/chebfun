function pass = test_ballfun()

% Example 1 :
S = [21,53,17];
f = ballfun(@(x,y,z)x,'cart',S);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
pass(1) = isequal(f,exact);

% Example 2 :
S = [42,12,23];
f = ballfun(@(x,y,z)y,'cart',S);
exact = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),S);
pass(2) = isequal(f,exact);

% Example 3 :
S = [20,31,42];
f = ballfun(@(x,y,z)z,'cart',S);
exact = ballfun(@(r,lam,th)r.*cos(th),S);
pass(3) = isequal(f,exact);

S = [50,50,50];

% Example 4 :
f = ballfun(@(x,y,z)x.*z,'cart',S);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*cos(th),S);
pass(4) = isequal(f,exact);

% Example 5 :
f = ballfun(@(x,y,z)sin(x.*y.*z),'cart',S);
exact = ballfun(@(r,lam,th)sin(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)),S);
pass(5) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end
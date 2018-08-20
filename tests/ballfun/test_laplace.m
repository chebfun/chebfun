function pass = test_laplace( )

% Test with different parity of m,p
% Example 1:
S = [38,37,40];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplace(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(1) = isequal(g,exact);

% Example 2:
S = [39,37,40];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplace(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(2) = isequal(g,exact);

% Example 3:
S = [38,37,41];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplace(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(3) = isequal(g,exact);

% Example 4:
S = [39,38,41];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplace(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(4) = isequal(g,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

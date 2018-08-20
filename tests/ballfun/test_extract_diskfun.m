function pass = test_extract_diskfun( )

eps = 1e-10;
S = [41,42,43];
m = S(1); n = S(2);

% Example 1
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(lam).*sin(th).^2,S);
g = extract_diskfun(f,'x','y');
h = diskfun(@(x,y)x.*y,[m,n]);
error = coeffs2(g-h,n,m);
pass(1) = max(abs(error(:))) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

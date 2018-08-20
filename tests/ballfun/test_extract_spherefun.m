function pass = test_extract_spherefun( )

eps = 1e-10;
S = [20,21,22];
n = S(2); p = S(3);

% Example 1
f = ballfun(@(r,lam,th)cos(lam).*sin(th),S);
g = extract_spherefun(f);
h = spherefun(@(lam,th)cos(lam).*sin(th));
error = coeffs2(g-h,p,n);
pass(1) = max(abs(error(:))) < eps;

% Example 2
f = ballfun(@(r,lam,th)r,[20,21,22]);
g = extract_spherefun(f, 0.5);
h = spherefun(@(lam,th)0.5,p,n);
error = coeffs2(g-h,p,n);
pass(2) = max(abs(error(:))) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

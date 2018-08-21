function pass = test_sum( ) 
eps = 1e-10;
S = [20,21,22];
n = S(2); p = S(3);

% Example 1
f = ballfun(@(r,lam,th)exp(r).*cos(lam).*sin(th),S);
g = sum(f);
exact = spherefun(@(lam,th)(exp(1)-1)*cos(lam).*sin(th),'vectorize');
error = coeffs2(g-exact,p,n);
pass(1) = max(abs(error(:))) < eps;

% Example 2
f = ballfun(@(r,lam,th)1,S);
g = sum(f);
exact = spherefun(@(lam,th)1,'vectorize');
error = coeffs2(g-exact,p,n);
pass(2) = max(abs(error(:))) < eps;

% Example 3 
f = ballfun(@(r,lam,th)cos(r),S);
g = sum(f);
exact = spherefun(@(lam,th)sin(1),'vectorize');
error = coeffs2(g-exact,p,n);
pass(3) = max(abs(error(:))) < eps;

% Example 4 
f = ballfun(@(r,lam,th)cos(r)+sin(r),S);
g = sum(f);
exact = spherefun(@(lam,th)sin(1)-cos(1)+1,'vectorize');
error = coeffs2(g-exact,p,n);
pass(4) = max(abs(error(:))) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

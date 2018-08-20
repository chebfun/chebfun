function pass = test_feval( )

S = [38,37,40];
zero = cheb.galleryballfun('zero',S);
r = 0.5; lam = pi; th = pi/2;
eps = 1e-15;

% Example 1:
f = ballfun(@(r,lam,th)1,S);
F = ballfunv(f,zero,zero);
vals = feval(F,r,lam,th);
vals_exact = zeros(1,1,1,3);
vals_exact(:,:,:,1) = 1;
pass(1) = norm(vals(:) - vals_exact(:),inf) < eps;

% Example 2:
f = ballfun(@(r,lam,th)r.*sin(th),S);
F = ballfunv(zero,f,zero);
vals = feval(F,r,lam,th);
vals_exact = zeros(1,1,1,3);
vals_exact(:,:,:,2) = 0.5;
pass(2) = norm(vals(:) - vals_exact(:),inf) < eps;

% Example 3:
f = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),S);
F = ballfunv(zero,zero,f);
vals = feval(F,r,lam,th);
vals_exact = zeros(1,1,1,3);
vals_exact(:,:,:,3) = -0.5;
pass(3) = norm(vals(:) - vals_exact(:),inf) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

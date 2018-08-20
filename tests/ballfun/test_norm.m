function pass = test_norm()

S = [30,30,30];
eps = 1e-10;

% Example 1
f = ballfun(@(r,lam,th)1,S);
norm_f = norm(f);
norm_exact = sqrt(4*pi/3);
pass(1) = abs(norm_f-norm_exact)<eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

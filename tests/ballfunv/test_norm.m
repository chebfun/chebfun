function pass = test_norm()

S = [30,30,30];
eps = 1e-10;

% Example 1
f = ballfun(@(r,lam,th)1,S);
v = ballfunv(f,f,f);
norm_v = norm(v);
norm_exact = sqrt(3*(4*pi/3));
pass(1) = abs(norm_v-norm_exact)<eps;

if (nargout > 0)
    pass = all(pass(:));
end
end

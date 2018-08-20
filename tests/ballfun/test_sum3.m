function pass = test_sum3( )
eps = 1e-10;

% Example 1
f = ballfun(@(r,lam,th)1,[20,20,20]);
I = sum3(f);
exact = 4*pi/3;
pass(1) = abs(I-exact)<eps;

% Example 2
f = ballfun(@(r,lam,th)cos(lam),[20,20,20]);
I = sum3(f);
exact = 0;
pass(2) = abs(I-exact)<eps;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*sin(lam),[20,20,20]);
I = sum3(f);
exact = 0;
pass(3) = abs(I-exact)<eps;

% Example 4
f = ballfun(@(r,lam,th)cos(th),[20,20,20]);
I = sum3(f);
exact = 0;
pass(4) = abs(I-exact)<eps;

% Example 5
f = ballfun(@(r,lam,th)sin(th),[20,20,20]);
I = sum3(f);
exact = pi^2/3;
pass(5) = abs(I-exact)<eps;

% Example 6
f = ballfun(@(r,lam,th)r.*sin(th).*(1+sin(lam)),[20,20,20]);
I = sum3(f);
exact = pi^2/4;
pass(6) = abs(I-exact)<eps;

if (nargout > 0)
    pass = all(pass(:));
end

end

function pass = test_sum2( )

S = [20,21,22];
m = S(1);

% Example 1
f = ballfun(@(r,lam,th)1,S);
ff = sum2(f);
g = chebfun(@(r)4*pi*r.^2,m);
pass(1) = max(abs(chebcoeffs(ff-g,m)))<1e-9;

% Example 2
f = ballfun(@(r,lam,th)cos(lam),S);
ff = sum2(f);
g = chebfun(@(r)0,m);
pass(2) = max(abs(chebcoeffs(ff-g,m)))<1e-9;

% Example 3
f = ballfun(@(r,lam,th)r.^2.*sin(lam),S);
ff = sum2(f);
g = chebfun(@(r)0,m);
pass(3) = max(abs(chebcoeffs(ff-g,m)))<1e-9;

if (nargout > 0)
    pass = all(pass(:));
end

end

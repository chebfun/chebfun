function pass = test_cosh( ) 
S = [20,21,22];

% Example 1
f = cosh(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)cosh(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

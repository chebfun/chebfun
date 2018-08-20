function pass = test_exp( ) 
S = [20,21,22];

% Example 1
f = exp(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)exp(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

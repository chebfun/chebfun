function pass = test_cos( ) 
S = [20,21,22];

% Example 1
f = cos(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)cos(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

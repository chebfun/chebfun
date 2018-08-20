function pass = test_sqrt( ) 
S = [20,21,22];

% Example 1
f = sqrt(ballfun(@(r,lam,th)r.^2,S));
exact = ballfun(@(r,lam,th)abs(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

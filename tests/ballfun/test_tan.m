function pass = test_tan( ) 
S = [20,21,22];

% Example 1
f = tan(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)tan(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

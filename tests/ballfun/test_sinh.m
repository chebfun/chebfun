function pass = test_sinh( ) 
S = [20,21,22];

% Example 1
f = sinh(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)sinh(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_sin( ) 
S = [20,21,22];

% Example 1
f = sin(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)sin(r),S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

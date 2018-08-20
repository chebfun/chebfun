function pass = test_log( ) 
S = [20,21,22];

% Example 1
f = log(ballfun(@(r,lam,th)exp(1i*th),S));
exact = ballfun(@(r,lam,th)1i*th,S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

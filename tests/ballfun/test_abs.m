function pass = test_abs( ) 
S = [20,21,22];

% Example 1
f = abs(ballfun(@(r,lam,th)exp(1i*lam),S));
exact = ballfun(@(r,lam,th)1+0*r,S);
pass(1) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

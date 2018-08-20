function pass = test_imag( ) 
S = [20,21,22];

% Example 1
f = imag(ballfun(@(r,lam,th)exp(1i*lam),S));
exact = ballfun(@(r,lam,th)sin(lam),S);
pass(1) = isequal(f,exact);

% Example 2
f = imag(ballfun(@(r,lam,th)exp(1i*th.*lam),S));
exact = ballfun(@(r,lam,th)sin(lam.*th),S);
pass(2) = isequal(f,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

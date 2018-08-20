function pass = test_minus( ) 
S = [18,19,20];

% Example 1
f = cheb.galleryballfun('random',S);
g = 2*f;
pass(1) = isequal(g-f,f);


% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(th),S);
g = f-5;
exact = ballfun(@(r,lam,th)r.^2.*cos(th)-5,S);
pass(2) = isequal(g,exact);

if (nargout > 0)
    pass = all(pass(:));
end
end

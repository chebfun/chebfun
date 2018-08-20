function pass = test_times( )

S = [38,37,40];
f = ballfun(@(r,lam,th)r.*cos(lam),S);
zero = cheb.galleryballfun('zero',S);

% Example 1:
F = ballfunv(f,zero,zero);
G = ballfunv(zero,f,zero);
H = times(F,G);
Hexact = ballfunv(zero,zero,zero);
pass(1) = isequal(H,Hexact);

% Example 2:
F = ballfunv(f,zero,zero);
G = ballfunv(-f,zero,zero);
H = times(F,G);
Hexact = ballfunv(-power(f,2),zero,zero);
pass(2) = isequal(H,Hexact);

% Example 3:
F = ballfunv(f,2*f,3*f);
G = ballfunv(-2*f,f,-f);
H = times(F,G);
Hexact = ballfunv(-2*power(f,2),2*power(f,2),-3*power(f,2));
pass(3) = isequal(H,Hexact);

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_cross( )

S = [38,37,40];
f = ballfun(@(r,lam,th)1,S);
zero = cheb.galleryballfun('zero',S);

% Example 1:
F = ballfunv(f,zero,zero);
G = ballfunv(zero,f,zero);
H = cross(F,G);
Hexact = ballfunv(zero,zero,f);
pass(1) = isequal(H,Hexact);

% Example 2:
F = ballfunv(f,zero,zero);
G = ballfunv(f,zero,zero);
H = cross(F,G);
Hexact = ballfunv(zero,zero,zero);
pass(2) = isequal(H,Hexact);

% Example 3:
F = ballfunv(zero,zero,f);
G = ballfunv(zero,f,zero);
H = cross(F,G);
Hexact = ballfunv(-f,zero,zero);
pass(3) = isequal(H,Hexact);

if (nargout > 0)
    pass = all(pass(:));
end
end

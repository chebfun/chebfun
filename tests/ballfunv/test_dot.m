function pass = test_dot( )

S = [38,37,40];
f = ballfun(@(r,lam,th)1,S);
zero = cheb.galleryballfun('zero',S);

% Example 1:
F = ballfunv(f,zero,zero);
G = ballfunv(zero,f,zero);
H = dot(F,G);
Hexact = zero;
pass(1) = isequal(H,Hexact);

% Example 2:
F = ballfunv(f,zero,zero);
G = ballfunv(-f,zero,zero);
H = dot(F,G);
Hexact = -power(f,2);
pass(2) = isequal(H,Hexact);

% Example 3:
F = ballfunv(f,2*f,3*f);
G = ballfunv(-2*f,f,-f);
H = dot(F,G);
Hexact = -3*power(f,2);
pass(3) = isequal(H,Hexact);

if (nargout > 0)
    pass = all(pass(:));
end
end

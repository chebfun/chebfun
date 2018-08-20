function pass = test_power( )

S = [38,37,40];
f = ballfun(@(r,lam,th)r.*cos(lam),S);
zero = cheb.galleryballfun('zero',S);
one = ballfun(@(r,lam,th)1,S);

% Example 1:
F = ballfunv(f,zero,zero);
n = 2;
H = power(F,n);
Hexact = ballfunv(power(f,n),zero,zero);
pass(1) = isequal(H,Hexact);

% Example 3:
F = ballfunv(f,2*f,3*f);
n = 0;
H = power(F,n);
Hexact = ballfunv(one,one,one);
pass(2) = isequal(H,Hexact);

% Example 3:
F = ballfunv(f,2*f,3*f);
n = 3;
H = power(F,n);
Hexact = ballfunv(power(f,3),8*power(f,3),27*power(f,3));
pass(3) = isequal(H,Hexact);

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_times( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

S = [38,37,40];
f = ballfun(@(r,lam,th)r.*cos(lam),S);
Z = cheb.galleryballfun('zero',S);

% Example 1:
F = ballfunv(f,Z,Z);
G = ballfunv(Z,f,Z);
H = times(F,G);
Hexact = ballfunv(Z,Z,Z);
pass(1) = norm(H - Hexact) < tol;

% Example 2:
F = ballfunv(f,Z,Z);
G = ballfunv(-f,Z,Z);
H = times(F,G);
Hexact = ballfunv(-power(f,2),Z,Z);
pass(2) = norm(H - Hexact) < tol;

% Example 3:
F = ballfunv(f,2*f,3*f);
G = ballfunv(-2*f,f,-f);
H = times(F,G);
Hexact = ballfunv(-2*power(f,2),2*power(f,2),-3*power(f,2));
pass(3) = norm(H - Hexact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

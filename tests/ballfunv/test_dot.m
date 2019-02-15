function pass = test_dot( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

f = ballfun(@(r,lam,th)1, 'spherical');
zero = ballfun(@(x,y,z)0);

% Example 1:
F = ballfunv(f,zero,zero);
G = ballfunv(zero,f,zero);
H = dot(F,G);
Hexact = zero;
pass(1) = norm(H-Hexact)<tol;

% Example 2:
F = ballfunv(f,zero,zero);
G = ballfunv(-f,zero,zero);
H = dot(F,G);
Hexact = -power(f,2);
pass(2) = norm(H-Hexact)<tol;

% Example 3:
F = ballfunv(f,2*f,3*f);
G = ballfunv(-2*f,f,-f);
H = dot(F,G);
Hexact = -3*power(f,2);
pass(3) = norm(H-Hexact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

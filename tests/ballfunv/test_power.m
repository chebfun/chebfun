function pass = test_power( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

f = ballfun(@(x,y,z)cos(x.*y).*z);
zero = ballfun(@(x,y,z)0);
one = ballfun(@(r,lam,th)1);

% Example 1:
F = ballfunv(f,zero,zero);
n = 2;
H = power(F,n);
Hexact = ballfunv(power(f,n),zero,zero);
pass(1) = norm(H-Hexact) < tol;

% Example 3:
F = ballfunv(f,2*f,3*f);
n = 0;
H = power(F,n);
Hexact = ballfunv(one,one,one);
pass(2) = norm(H-Hexact) < tol;

% Example 3:
F = ballfunv(f,2*f,3*f);
n = 3;
H = power(F,n);
Hexact = ballfunv(power(f,3),8*power(f,3),27*power(f,3));
pass(3) = norm(H-Hexact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

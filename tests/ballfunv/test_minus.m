function pass = test_minus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f1 = ballfun(@(x,y,z)sin(y));
f2 = ballfun(@(x,y,z)cos(x));
f3 = ballfun(@(x,y,z)sin(y.*z));
f = ballfunv(f1,f2,f3);
g = 2*f;
pass(1) = norm(g-f - f)<tol;

% Example 2
f1 = ballfun(@(x,y,z)y);
f2 = ballfun(@(x,y,z)cos(x));
f3 = ballfun(@(x,y,z)z.*x);
f = ballfunv(f1,f2,f3);
g = f-5;
e1 = ballfun(@(x,y,z)y-5);
e2 = ballfun(@(x,y,z)cos(x)-5);
e3 = ballfun(@(x,y,z)z.*x-5);
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

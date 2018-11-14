function pass = test_plus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f1 = ballfun(@(x,y,z)cos(y));
f2 = ballfun(@(x,y,z)sin(x));
f3 = ballfun(@(x,y,z)cos(y.*z));
f = ballfunv(f1,f2,f3);
g = f+f;
e1 = ballfun(@(x,y,z)2*cos(y));
e2 = ballfun(@(x,y,z)2*sin(x));
e3 = ballfun(@(x,y,z)2*cos(y.*z));
exact = ballfunv(e1,e2,e3);
pass(1) = norm(g-exact)<tol;

% Example 2
f1 = ballfun(@(x,y,z)sin(y));
f2 = ballfun(@(x,y,z)cos(x));
f3 = ballfun(@(x,y,z)sin(y.*z));
f = ballfunv(f1,f2,f3);
g = f+1;
e1 = ballfun(@(x,y,z)sin(y)+1);
e2 = ballfun(@(x,y,z)cos(x)+1);
e3 = ballfun(@(x,y,z)sin(y.*z)+1);
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

% Example 3
f1 = ballfun(@(x,y,z)y);
f2 = ballfun(@(x,y,z)cos(x));
f3 = ballfun(@(x,y,z)z.*x);
f = ballfunv(f1,f2,f3);
g = 3+f;
e1 = ballfun(@(x,y,z)y+3);
e2 = ballfun(@(x,y,z)cos(x)+3);
e3 = ballfun(@(x,y,z)z.*x+3);
exact = ballfunv(e1,e2,e3);
pass(3) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

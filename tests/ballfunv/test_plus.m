function pass = test_plus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f1 = ballfun(@(x,y,z)cos(y),'cart');
f2 = ballfun(@(x,y,z)sin(x),'cart');
f3 = ballfun(@(x,y,z)cos(y.*z),'cart');
f = ballfunv(f1,f2,f3);
g = f+f;
e1 = ballfun(@(x,y,z)2*cos(y),'cart');
e2 = ballfun(@(x,y,z)2*sin(x),'cart');
e3 = ballfun(@(x,y,z)2*cos(y.*z),'cart');
exact = ballfunv(e1,e2,e3);
pass(1) = norm(g-exact)<tol;

% Example 2
f1 = ballfun(@(x,y,z)sin(y),'cart');
f2 = ballfun(@(x,y,z)cos(x),'cart');
f3 = ballfun(@(x,y,z)sin(y.*z),'cart');
f = ballfunv(f1,f2,f3);
g = f+1;
e1 = ballfun(@(x,y,z)sin(y)+1,'cart');
e2 = ballfun(@(x,y,z)cos(x)+1,'cart');
e3 = ballfun(@(x,y,z)sin(y.*z)+1,'cart');
exact = ballfunv(e1,e2,e3);
pass(2) = norm(g-exact)<tol;

% Example 3
f1 = ballfun(@(x,y,z)y,'cart');
f2 = ballfun(@(x,y,z)cos(x),'cart');
f3 = ballfun(@(x,y,z)z.*x,'cart');
f = ballfunv(f1,f2,f3);
g = 3+f;
e1 = ballfun(@(x,y,z)y+3,'cart');
e2 = ballfun(@(x,y,z)cos(x)+3,'cart');
e3 = ballfun(@(x,y,z)z.*x+3,'cart');
exact = ballfunv(e1,e2,e3);
pass(3) = norm(g-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

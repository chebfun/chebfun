function pass = test_real( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)z,'cart');
V = real(ballfunv(f,f,f));
g = ballfun(@(x,y,z)z,'cart');
exact = ballfunv(g,g,g);
pass(1) = norm(V-exact)<tol;

% Example 2
f = ballfun(@(x,y,z)y+1i*z,'cart');
V = real(ballfunv(f,f,f));
g = ballfun(@(x,y,z)y,'cart');
exact = ballfunv(g,g,g);
pass(2) = norm(V-exact)<tol;

% Example 3
f1 = ballfun(@(x,y,z)x,'cart');
f2 = ballfun(@(x,y,z)1i*z,'cart');
f3 = ballfun(@(x,y,z)cos(y)+1i*sin(x),'cart');
V = real(ballfunv(f1,f2,f3));
g1 = ballfun(@(x,y,z)x,'cart');
g2 = ballfun(@(x,y,z)0,'cart');
g3 = ballfun(@(x,y,z)cos(y),'cart');
exact = ballfunv(g1,g2,g3);
pass(3) = norm(V-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

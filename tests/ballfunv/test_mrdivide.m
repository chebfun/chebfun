function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

zero = ballfun(@(x,y,z)0);

% Example 1
f = ballfun(@(x,y,z)cos(y.*z));
V = ballfunv(f,zero,zero);
W = V/2;
g = ballfun(@(x,y,z)cos(y.*z)/2);
exact = ballfunv(g,zero,zero);
pass(1) = norm(W-exact)<tol;

% Example 2
f = ballfun(@(x,y,z)sin(y.*z));
V = ballfunv(f,f,f);
W = V/(-3);
g = ballfun(@(x,y,z)-sin(y.*z)/3);
exact = ballfunv(g,g,g);
pass(2) = norm(W-exact)<tol;

% Example 3
f = ballfun(@(x,y,z)2*sin(x));
V = ballfunv(f,zero,f);
W = V/1i;
g = ballfun(@(x,y,z)2*sin(x)/1i);
exact = ballfunv(g,zero,g);
pass(3) = norm(W-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

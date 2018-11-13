function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)x+y,'cart');
g = f/2;
exact = ballfun(@(x,y,z)(x+y)/2,'cart');
pass(1) = norm( g - exact ) < tol;

% Example 2
f = ballfun(@(x,y,z)x.*z,'cart');
g = f/(-3);
exact = ballfun(@(x,y,z)-x.*z/3,'cart');
pass(2) = norm( g - exact ) < tol;

% Example 3
f = ballfun(@(x,y,z)y,'cart');
g = f/1i;
exact = ballfun(@(x,y,z)-1i*y,'cart');
pass(3) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

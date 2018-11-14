function pass = test_mrdivide( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)x+y);
g = f/2;
exact = ballfun(@(x,y,z)(x+y)/2);
pass(1) = norm( g - exact ) < tol;

% Example 2
f = ballfun(@(x,y,z)x.*z);
g = f/(-3);
exact = ballfun(@(x,y,z)-x.*z/3);
pass(2) = norm( g - exact ) < tol;

% Example 3
f = ballfun(@(x,y,z)y);
g = f/1i;
exact = ballfun(@(x,y,z)-1i*y);
pass(3) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

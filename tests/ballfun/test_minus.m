function pass = test_minus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)x+y+x.^2);
g = 2*f;
h = g - f; 
pass(1) = norm( h - f ) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.*cos(th), 'spherical');
g = f-5;
exact = ballfun(@(r,lam,th)r.*cos(th)-5, 'spherical');
pass(2) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

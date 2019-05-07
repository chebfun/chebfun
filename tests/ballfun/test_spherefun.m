function pass = test_spherefun( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)cos(x.*y));
g = spherefun(f);
h = spherefun(@(x,y,z)cos(x.*y));
pass(1) = norm(g-h) < tol;

% Example 2
f = ballfun(@(x,y,z)sin(z));
g = spherefun(f, 0.5);
h = spherefun(@(x,y,z)sin(0.5*z));
pass(2) = norm(g-h) < tol;

% Example 3
f = ballfun(@(x,y,z)x);
g = f(10,:,:,'spherical');
h = spherefun(@(x,y,z)10*x);
pass(3) = norm(g-h) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

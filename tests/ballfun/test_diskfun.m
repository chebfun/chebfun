function pass = test_diskfun( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)cos(x.*y));
g = diskfun(f);
h = diskfun(@(x,y)cos(x.*y));
pass(1) = norm(g-h) < tol;

% Example 2
f = ballfun(@(x,y,z)cos(x.*y));
g = diskfun(f,'z');
h = diskfun(@(x,y)cos(x.*y));
pass(2) = norm(g-h) < tol;

% Example 3
f = ballfun(@(x,y,z)x.*z);
g = diskfun(f,'y');
h = diskfun(@(x,y)x.*y);
pass(3) = norm(g-h) < tol;

% Example 4
f = ballfun(@(x,y,z)x.*y);
g = diskfun(f,'z');
h = diskfun(@(x,y)x.*y);
pass(4) = norm(g-h) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

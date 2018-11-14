function pass = test_tan( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = tan(ballfun(@(x,y,z)sin(y)));
exact = ballfun(@(x,y,z)tan(sin(y)));
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

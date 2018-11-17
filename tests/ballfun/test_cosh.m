function pass = test_cosh( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = cosh(ballfun(@(x,y,z)y));
exact = ballfun(@(x,y,z)cosh(y));
pass(1) = norm( f - exact) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

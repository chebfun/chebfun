function pass = test_sin( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = sin(ballfun(@(r,lam,th)r));
exact = ballfun(@(r,lam,th)sin(r));
pass(1) = norm( f - exact ) < tol; 

if (nargout > 0)
    pass = all(pass(:));
end
end

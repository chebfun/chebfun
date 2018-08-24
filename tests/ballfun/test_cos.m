function pass = test_cos( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Example 1
f = cos(ballfun(@(r,lam,th)r));
exact = ballfun(@(r,lam,th)cos(r));
pass(1) = norm( f - exact ) < tol ;

if (nargout > 0)
    pass = all(pass(:));
end
end

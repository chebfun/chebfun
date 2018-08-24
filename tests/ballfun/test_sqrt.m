function pass = test_sqrt( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps; 

% Example 1
f = sqrt(ballfun(@(r,lam,th)r.^2));
exact = ballfun(@(r,lam,th)abs(r));
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_conj( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = conj(ballfun(@(x,y,z)x+1i*y));
exact = ballfun(@(x,y,z)x-1i*y);
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

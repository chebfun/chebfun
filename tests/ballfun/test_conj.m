function pass = test_conj( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = conj(ballfun(@(r,lam,th)exp(1i*lam)));
exact = ballfun(@(r,lam,th)exp(-1i*lam));
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

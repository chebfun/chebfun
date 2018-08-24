function pass = test_abs( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Example 1
f = abs(ballfun(@(r,lam,th)exp(1i*lam)));
exact = ballfun(@(r,lam,th)1+0*r);
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

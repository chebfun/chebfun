function pass = test_real( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [20,21,22];

% Example 1
f = real(ballfun(@(r,lam,th)exp(1i*lam),S));
exact = ballfun(@(r,lam,th)cos(lam),S);
pass(1) = norm( f - exact ) < tol;

% Example 2
f = real(ballfun(@(r,lam,th)exp(1i*th.*lam),S));
exact = ballfun(@(r,lam,th)cos(lam.*th),S);
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

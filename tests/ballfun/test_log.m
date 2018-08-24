function pass = test_log( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [20,21,22];

% Example 1
f = log(ballfun(@(r,lam,th)exp(1i*th),S));
exact = ballfun(@(r,lam,th)1i*th,S);
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

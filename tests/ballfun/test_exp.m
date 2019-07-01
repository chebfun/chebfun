function pass = test_exp( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = exp(ballfun(@(x,y,z)z));
exact = ballfun(@(x,y,z)exp(z));
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

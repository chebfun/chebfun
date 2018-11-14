function pass = test_sqrt( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e8*pref.techPrefs.chebfuneps; 

% Example 1
f = sqrt(ballfun(@(x,y,z)2));
exact = ballfun(@(x,y,z)sqrt(2));
pass(1) = norm( f - exact ) < tol;

% Example 1
f = sqrt(ballfun(@(x,y,z)x.^4));
exact = ballfun(@(x,y,z)x.^2);
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

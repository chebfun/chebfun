function pass = test_cos( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Example 1
f = cos(ballfun(@(x,y,z)x,'cart'));
exact = ballfun(@(x,y,z)cos(x),'cart');
pass(1) = norm( f - exact ) < tol ;

if (nargout > 0)
    pass = all(pass(:));
end
end

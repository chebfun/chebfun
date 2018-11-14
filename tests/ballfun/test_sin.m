function pass = test_sin( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = sin(ballfun(@(x,y,z)x.*y.*z));
exact = ballfun(@(x,y,z)sin(x.*y.*z));
pass(1) = norm( f - exact ) < tol; 

if (nargout > 0)
    pass = all(pass(:));
end
end

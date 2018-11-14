function pass = test_tanh( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps; 

% Example 1
f = tanh(ballfun(@(x,y,z)cos(x.*z)));
exact = ballfun(@(x,y,z)tanh(cos(x.*z)));
pass(1) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

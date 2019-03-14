function pass = test_isempty( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun();
pass(1) = isempty(f);

% Example 2
f = ballfun([]);
pass(2) = isempty(f);

if (nargout > 0)
    pass = all(pass(:));
end
end
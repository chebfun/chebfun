function pass = test_size( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)1);
b1 = size(f) == [1,1,1];
pass(1) = all(b1);

% Example 2
f = ballfun(@(x,y,z)x);
b1 = size(f) == [2,3,3];
pass(2) = all(b1);

if (nargout > 0)
    pass = all(pass(:));
end
end

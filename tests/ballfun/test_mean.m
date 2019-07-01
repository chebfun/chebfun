function pass = test_mean( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = pref.techPrefs.chebfuneps;

% Example 1:
f = ballfun(@(x,y,z)1);
I = mean(f,1);
g = spherefun(1);
pass(1) = ( norm(I-g) < tol );

% Example 2:
f = ballfun(@(x,y,z)2);
I = mean(f,2);
g = diskfun(2);
pass(2) = ( norm(I-g) < tol );

% Example 3:
f = ballfun(@(x,y,z)3);
I = mean(f,3);
g = diskfun(3);
pass(3) = ( norm(I-g) < tol );
end
function pass = test_mean2( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1:
f = ballfun(@(x,y,z)1);
I = mean2(f,[1,2]);
g = chebfun(@(x)1,[-pi,pi],'trig');
pass(1) = ( norm(I-g) < tol );

% Example 2:
f = ballfun(@(x,y,z)2);
I = mean2(f,[2,3]);
g = chebfun(@(x)2);
pass(2) = ( norm(I-g) < tol );

% Example 3:
f = ballfun(@(x,y,z)3);
I = mean2(f,[1,3]);
g = chebfun(@(x)3,[-pi,pi],'trig');
pass(3) = ( norm(I-g) < tol );
end
function pass = test_sample( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(x,y,z)x);
X = sample(f);
Y = zeros(2,3,3);
Y(2,:,2) = [-1,0.5,0.5];
pass(1) = norm( X(:)-Y(:) ) < tol;

% Example 2
f = ballfun(@(x,y,z)z);
X = sample(f, 3, 1, 4);
r = [0;sqrt(0.5);1];
cosT = [1,0.5,-0.5,-1];
Y = reshape(r*cosT,3,1,4);
pass(2) = norm( X(:)-Y(:) ) < tol;

% Example 3
f = ballfun(@(x,y,z)z);
pass(3) = norm(ballfun(sample(f,4,4,4))-f) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

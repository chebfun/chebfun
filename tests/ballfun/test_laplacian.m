function pass = test_laplacian( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e6*pref.techPrefs.chebfuneps;

% Test with different parity of m,p
% Example 1:
f = ballfun(@(x,y,z)x.^2+y.^2+z.^2);
g = laplacian(f);
exact = ballfun(@(x,y,z)6);
pass(1) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_abs( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Example 1
f = abs(ballfun(@(x,y,z)1i*x,'cart'));
exact = ballfun(@(x,y,z)abs(x),'cart');
pass(1) = norm( f - exact ) < tol;

% Example 2
f = abs(ballfun(@(x,y,z)x+1i*y,'cart'));
exact = ballfun(@(x,y,z)sqrt(x.^2+y.^2),'cart');
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

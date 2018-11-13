function pass = test_real( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = real(ballfun(@(x,y,z)x+1i*y,'cart'));
exact = ballfun(@(x,y,z)x,'cart');
pass(1) = norm( f - exact ) < tol;

% Example 2
f = real(ballfun(@(x,y,z)sin(z)+1i*cos(y),'cart'));
exact = ballfun(@(x,y,z)sin(z),'cart');
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

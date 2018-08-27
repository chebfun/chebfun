function pass = test_log( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = log(ballfun(@(r,lam,th)10));
exact = ballfun(@(r,lam,th)log(10));
pass(1) = norm( f - exact ) < tol;

% Example 2
f = log(ballfun(@(r,lam,th)exp(r)));
exact = ballfun(@(r,lam,th)r);
pass(2) = norm( f - exact ) < tol;

% Example 3
f = log(ballfun(@(x,y,z)exp(x.*z),'cart'));
exact = ballfun(@(x,y,z)x.*z,'cart');
pass(3) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

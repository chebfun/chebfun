function pass = test_log( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = log(ballfun(@(r,lam,th)10,'spherical'));
exact = ballfun(@(r,lam,th)log(10),'spherical');
pass(1) = norm( f - exact ) < tol;

% Example 2
f = log(ballfun(@(x,y,z)exp(y)));
exact = ballfun(@(x,y,z)y);
pass(2) = norm( f - exact ) < tol;

% Example 3
f = log(ballfun(@(x,y,z)exp(x.*z)));
exact = ballfun(@(x,y,z)x.*z);
pass(3) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

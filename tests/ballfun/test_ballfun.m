function pass = test_ballfun( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = pref.techPrefs.chebfuneps; 

% Example 1 :
f = ballfun(@(x,y,z)x, 'cart');
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam));
pass(1) = norm( f - exact ) < tol;

% Example 2 :
f = ballfun(@(x,y,z)y,'cart');
exact = ballfun(@(r,lam,th)r.*sin(th).*sin(lam));
pass(2) = norm( f - exact ) < tol;

% Example 3 :
f = ballfun(@(x,y,z)z,'cart');
exact = ballfun(@(r,lam,th)r.*cos(th));
pass(3) = norm( f - exact ) < tol;

% Example 4 :
f = ballfun(@(x,y,z)x.*z,'cart');
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*cos(th));
pass(4) = norm( f - exact ) < tol;

% Example 5 :
f = ballfun(@(x,y,z)sin(x.*y.*z),'cart');
exact = ballfun(@(r,lam,th)sin(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)));
pass(5) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
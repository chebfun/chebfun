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

% Example 6 : (Test 'coeffs' flag.) with rcos(th)
m = 10; n = 11; p = 12;
F = zeros(m,n,p);
F(2,floor(n/2)+1,floor(p/2))=1/2;F(2,floor(n/2)+1,floor(p/2)+2)=1/2;
f = ballfun(F, 'coeffs');
g = ballfun(@(r,lam,th)r.*cos(th));
pass(6) = norm( f - g ) < tol;

% Example 7 :
f = ballfun(@(x,y,z)cos(x.*y.*z),'cart');
g = ballfun(@(r,lam,th)cos(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)));
pass(7) = norm( f - g ) < tol;

% Example 8
f = ballfun(@(x,y,z)cos(x.*y.*z),'cart');
g = ballfun(@(x,y,z)cos(x.*y.*z),'cart','vectorize');
pass(8) = norm( f - g ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
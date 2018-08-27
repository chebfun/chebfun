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

% Example 6: (Test 'coeffs' flag.) 
m = 10; n = 11; p = 12;
F = zeros(m,n,p);
F(1,floor(n/2),floor(p/2))=1/4;F(1,floor(n/2),floor(p/2)+2)=1/4;
F(1,floor(n/2)+2,floor(p/2))=1/4;F(1,floor(n/2)+2,floor(p/2)+2)=1/4;
F(2,floor(n/2),floor(p/2))=1/4;F(2,floor(n/2),floor(p/2)+2)=1/4;
F(2,floor(n/2)+2,floor(p/2))=1/4;F(2,floor(n/2)+2,floor(p/2)+2)=1/4;

f = ballfun(F, 'coeffs');
g = ballfun(@(r,lam,th)(r+1).*cos(lam).*cos(th));

pass(6) = norm( f - g ) < 2*tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
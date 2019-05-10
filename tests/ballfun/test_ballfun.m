function pass = test_ballfun( pref )

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = pref.techPrefs.chebfuneps; 

% Example 1 :
f = ballfun(@(x,y,z)x);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam),'spherical');
pass(1) = norm( f - exact ) < tol;

% Example 2 :
f = ballfun(@(x,y,z)y);
exact = ballfun(@(r,lam,th)r.*sin(th).*sin(lam),'spherical');
pass(2) = norm( f - exact ) < tol;

% Example 3 :
f = ballfun(@(x,y,z)z);
exact = ballfun(@(r,lam,th)r.*cos(th),'spherical');
pass(3) = norm( f - exact ) < tol;

% Example 4 :
f = ballfun(@(x,y,z)x.*z);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*cos(th),'spherical');
pass(4) = norm( f - exact ) < tol;

% Example 5 :
f = ballfun(@(x,y,z)sin(x.*y.*z));
exact = ballfun(@(r,lam,th)sin(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)),'spherical');
pass(5) = norm( f - exact ) < tol;

% Example 6 : (Test 'coeffs' flag.) with rcos(th)
m = 10; n = 11; p = 12;
F = zeros(m,n,p);
F(2,floor(n/2)+1,floor(p/2))=1/2;F(2,floor(n/2)+1,floor(p/2)+2)=1/2;
f = ballfun(F, 'coeffs');
g = ballfun(@(r,lam,th)r.*cos(th),'spherical');
pass(6) = norm( f - g ) < tol;

% Example 7 :
f = ballfun(@(x,y,z)cos(x.*y.*z));
g = ballfun(@(r,lam,th)cos(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)),'spherical');
pass(7) = norm( f - g ) < tol;

% Example 8
f = ballfun(@(x,y,z)cos(x.*y.*z));
g = ballfun(@(x,y,z)cos(x.*y.*z),'vectorize');
pass(8) = norm( f - g ) < tol;

% Example 9 :
f = ballfun(@(x,y,z)x.*z);
exact = ballfun(@(r,lam,th)r.*sin(th).*cos(lam).*r.*cos(th),'spherical');
pass(9) = norm( f - exact ) < tol;

% Example 10 :
f = ballfun(@(x,y,z)sin(x.*y.*z));
exact = ballfun(@(r,lam,th)sin(r.*sin(th).*cos(lam).*r.*sin(th).*sin(lam).*r.*cos(th)),'spherical');
pass(10) = norm( f - exact ) < tol;

% Example 11:
f = ballfun(@(x,y,z)x*y,'vectorize');
exact = ballfun(@(x,y,z)x.*y);
pass(11) = norm( f - exact ) < tol;

% Example 12:
f = ballfun(@(x,y,z)x*z*y,'vectorize');
exact = ballfun(@(x,y,z)x.*z.*y);
pass(12) = norm( f - exact ) < tol;

% Example 13
f = ballfun(@(x,y,z)cos(x.*y),[50,51,52]);
exact = ballfun(@(x,y,z)cos(x.*y));
pass(13) = norm( f - exact ) < tol;

% Example 14
exact = ballfun(@(x,y,z)cos(x.*y));
f = ballfun(exact, [51,50,49], 'vectorize');
pass(14) = norm( f - exact ) < tol;

% Example 15
exact = ballfun(@(x,y,z)sin(x.*y));
E = coeffs3(exact,49,50,51);
f = ballfun(E, 'coeffs', [51,53,51]);
pass(15) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
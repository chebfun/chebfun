function pass = test_plus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(ones(21,20,22));
V1 = ballfun.coeffs2vals(f.coeffs);
g = f+f;
V2 = ballfun.coeffs2vals(g.coeffs);
pass(1) = ( norm(V2(:)-2*V1(:),inf) < tol );

% Example 2
f = ballfun(@(x,y,z)x.^2.*cos(y)-1);
g = f+1;
exact = ballfun(@(x,y,z)x.^2.*cos(y));
pass(2) = norm( g - exact ) < tol;

% Example 3
f = ballfun(@(x,y,z)y.*sin(z));
g = 3+f;
exact = ballfun(@(x,y,z)y.*sin(z)+3);
pass(3) = norm( g - exact ) < tol;

% Example 4
f = ballfun(@(x,y,z)x.*sin(z).^2.*cos(y));
g = 3+f+2;
exact = ballfun(@(x,y,z)x.*sin(z).^2.*cos(y)+5);
pass(4) = norm( g - exact ) < tol;

% Example 5
f = ballfun(@(r,lam,th)1);
g = f+f;
exact = ballfun(@(r,lam,th)2);
pass(5) = norm( g - exact ) < tol;

% Example 6
f = ballfun(0) + 1;
exact = ballfun(1);
pass(6) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end

end

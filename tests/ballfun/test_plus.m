function pass = test_plus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [19,20,21];
% Example 1
f = ballfun(ones(20,21,22));
V1 = ballfun.coeffs2vals(f.coeffs);
g = f+f;
V2 = ballfun.coeffs2vals(g.coeffs);
pass(1) = ( norm(V2(:)-2*V1(:),inf) < tol );

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(th)+1,S);
g = f+1;
exact = ballfun(@(r,lam,th)r.^2.*cos(th)+2,S);
pass(2) = norm( g - exact ) < tol;

% Example 3
f = ballfun(@(r,lam,th)r.*sin(th)+1,S);
g = 3+f;
exact = ballfun(@(r,lam,th)r.*sin(th)+4,S);
pass(3) = norm( g - exact ) < tol;

% Example 4
f = ballfun(@(r,lam,th)r.*sin(th).^2.*cos(lam)+1,S);
g = 3+f+2;
exact = ballfun(@(r,lam,th)r.*sin(th).^2.*cos(lam)+6,S);
pass(4) = norm( g - exact ) < tol;

% Example 5
f = ballfun(@(r,lam,th)1);
g = f+f;
exact = ballfun(@(r,lam,th)2);
pass(5) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end

end

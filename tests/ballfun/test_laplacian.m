function pass = test_laplacian( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with different parity of m,p
% Example 1:
S = [38,37,40];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplacian(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(1) = norm( g - exact ) < tol;

% Example 2:
S = [39,37,40];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplacian(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(2) = norm( g - exact ) < tol;

% Example 3:
S = [38,37,41];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplacian(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(3) = norm( g - exact ) < tol;

% Example 4:
S = [39,38,41];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th).^2,S);
g = laplacian(f);
exact = ballfun(@(r,lam,th)6*cos(lam).*sin(th).^2+cos(lam).*(4*cos(th).^2-2*sin(th).^2)+...
                                -cos(lam),S);
pass(4) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

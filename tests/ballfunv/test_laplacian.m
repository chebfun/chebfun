function pass = test_laplacian( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e7*pref.techPrefs.chebfuneps;

% Test with different parity of m,p
% Example 1:
f = ballfun(@(x,y,z)x.^2+y.^2+z.^2);
V = ballfunv(f,2*f,-f);
g = laplacian(V);
exact = ballfunv(ballfun(@(x,y,z)6),ballfun(@(x,y,z)12),ballfun(@(x,y,z)-6));
pass(1) = norm( g - exact ) < tol;

% Example 2:
f1 = ballfun(@(x,y,z)cos(x.*y));
f2 = ballfun(@(x,y,z)sin(y.*z));
f3 = ballfun(@(x,y,z)z.^3);
v = ballfunv(f1,f2,f3);
g = laplacian(v);
exact = grad(div(v))-curl(curl(v));
pass(2) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

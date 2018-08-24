function pass = test_minus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [18,19,20];

% Example 1
f = cheb.galleryballfun('random',S);
g = 2*f;
h = g - f; 
pass(1) = norm( h - f ) < tol;

% Example 2
f = ballfun(@(r,lam,th)r.^2.*cos(th),S);
g = f-5;
exact = ballfun(@(r,lam,th)r.^2.*cos(th)-5,S);
pass(2) = norm( g - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

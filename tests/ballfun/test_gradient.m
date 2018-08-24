function pass = test_gradient( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1: xyz
S = [38,37,40];
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th).*r.*cos(th),S);
v = gradient(f);
exactx = ballfun(@(r,lam,th)r.*sin(lam).*sin(th).*r.*cos(th),S);
exacty = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*cos(th),S);
exactz = ballfun(@(r,lam,th)r.*cos(lam).*sin(th).*r.*sin(lam).*sin(th),S);
exact = ballfunv(exactx,exacty,exactz);
pass(1) = norm( v - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end
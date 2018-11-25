function pass = test_norm( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)1, 'spherical');
v = ballfunv(f,f,f);
norm_v = norm(v);
norm_exact = sqrt(3*(4*pi/3));
pass(1) = abs(norm_v-norm_exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

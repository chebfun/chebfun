function pass = test_norm( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th) 1, 'spherical');
norm_f = norm( f );
norm_exact = sqrt(4*pi/3);
pass(1) = abs( norm_f-norm_exact )<tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

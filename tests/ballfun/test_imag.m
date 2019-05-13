function pass = test_imag( pref ) 
% Test imag

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = imag(ballfun(@(x,y,z)x+1i*y.*z));
exact = ballfun(@(x,y,z)y.*z);
pass(1) = norm( f - exact ) < tol;

% Example 2
f = imag(ballfun(@(x,y,z)y));
exact = ballfun(@(x,y,z)0);
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

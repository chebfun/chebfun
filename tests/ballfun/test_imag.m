function pass = test_imag( pref ) 
% Test the Helmholtz solver with Dirichlet BC

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = imag(ballfun(@(x,y,z)x+1i*y.*z,'cart'));
exact = ballfun(@(x,y,z)y.*z,'cart');
pass(1) = norm( f - exact ) < tol;

% Example 2
f = imag(ballfun(@(x,y,z)y,'cart'));
exact = ballfun(@(x,y,z)0,'cart');
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = test_imag( pref ) 
% Test the Helmholtz solver with Dirichlet BC

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [20,21,22];

% Example 1
f = imag(ballfun(@(r,lam,th)exp(1i*lam),S));
exact = ballfun(@(r,lam,th)sin(lam),S);
pass(1) = norm( f - exact ) < tol;

% Example 2
f = imag(ballfun(@(r,lam,th)exp(1i*th.*lam),S));
exact = ballfun(@(r,lam,th)sin(lam.*th),S);
pass(2) = norm( f - exact ) < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end

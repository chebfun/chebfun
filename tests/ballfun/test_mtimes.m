function pass = test_mtimes( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function coeffs 1
f = ballfun(ones(21,20,22));
g = ballfun(2*ones(21,20,22));

pass(1) = (norm(2*f-g)<tol && norm(f*2-g)<tol);
end

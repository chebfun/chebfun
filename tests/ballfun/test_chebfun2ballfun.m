function pass = test_chebfun2ballfun( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Test with function sin(r)
f = chebfun(@(r)sin(r));
g = ballfun.chebfun2ballfun(f);
h = ballfun(@(r,lam,th) sin(r));
pass(1) = norm( g - h ) < tol;

end

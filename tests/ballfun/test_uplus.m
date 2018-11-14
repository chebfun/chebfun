function pass = test_uplus( pref ) 
% Test with function rand = +rand

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

f = ballfun(@(x,y,z)cos(z).*y+sin(x));
g = +f;

pass(1) = norm( f - g ) < tol;
end

function pass = test_uplus( pref ) 
% Test with function rand = +rand

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

S = [22,23,14];
f = cheb.galleryballfun('random',S);
g = +f;

pass(1) = norm( f - g ) < tol;
end

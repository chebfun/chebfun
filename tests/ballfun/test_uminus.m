function pass = test_uminus( pref ) 
% Test with function rand : rand + -rand = 0

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

f = ballfun(@(x,y,z)cos(x).*z+sin(x));
g = -f;

pass(1) = norm(f+g) < tol;
end

function pass = test_size( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function coeffs 1
f = ballfun(ones(20,21,22));
g = ballfun(ones(1000,1,1));
b1 = size(f) == [20,21,22];
b2 = size(g) == [1000,1,1];

pass(1) = (min(b1)==1 && min(b2)==1);
end

function pass = test_isequal( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function 1
f = ballfun(ones(21,18,22));
F = ballfunv(f,f,f);
G = F+F-F;

pass(1) = (isequal(F,G) && isequal(G,F));
end

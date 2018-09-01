function pass = test_mtimes( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function coeffs 1
f = ballfun(ones(20,21,22));
F = ballfunv(f,f,f);
g = ballfun(2*ones(20,21,22));
G = ballfunv(g,g,g);
pass(1) = norm(2*F-G)<tol; 
pass(2) = norm(F*2-G)<tol;
end

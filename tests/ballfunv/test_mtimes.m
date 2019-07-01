function pass = test_mtimes( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Test with function coeffs 1
f = ballfun(ones(21,22,22));
F = ballfunv(f,f,f);
g = ballfun(2*ones(21,22,22));
G = ballfunv(g,g,g);
pass(1) = norm(2*F-G)<tol; 
pass(2) = norm(F*2-G)<tol;

% Example 2
% Multiply ballfunv by ballfun
V = ballfunv(@(x,y,z)x,@(x,y,z)y,@(x,y,z)z);
f = ballfun(@(x,y,z)cos(y));
exact = ballfunv(@(x,y,z)x.*cos(y),@(x,y,z)y.*cos(y),@(x,y,z)z.*cos(y));
pass(3) = norm(V.*f-exact) < tol;
pass(4) = norm(f.*V-exact) < tol;
end

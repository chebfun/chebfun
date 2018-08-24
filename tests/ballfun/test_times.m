function pass = test_times( pref ) 
% Test with function cos(cos(lam)*sin(th))

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function 1, 2*3=6
f = ballfun(2*ones(20,21,22),"vals");
g = ballfun(3*ones(20,21,22),"vals");
h = ballfun(6*ones(20,21,22),"vals");

pass(1) = ( norm(f*g-h)<tol );
pass(2) = ( norm(g*f-h)<tol ); 
end
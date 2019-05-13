function pass = test_times( pref ) 
% Test with function cos(cos(lam)*sin(th))

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function 1, 2*3=6
f = ballfun(2*ones(21,12,22));
g = ballfun(3*ones(23,20,24));
h = ballfun(6*ones(25,20,10));

pass(1) = ( norm(f*g-h)<tol );
pass(2) = ( norm(g*f-h)<tol ); 
end
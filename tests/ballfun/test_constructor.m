function pass = test_constructor( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Can we make a ballfun object:
n = 10; 
f = ballfun( ones(n,n,n) );
pass(1) = 1;

end

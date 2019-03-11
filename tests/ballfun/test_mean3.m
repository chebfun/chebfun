function pass = test_mean3( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function 1
f = ballfun(@(r,lam,th) 1, 'spherical');
I = mean3(f);
pass(1) = ( abs(I-1) < tol );
end
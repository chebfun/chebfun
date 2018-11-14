function pass = test_uminus( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

% Example 1
F1 = ballfun(@(x,y,z)cos(y));
F2 = ballfun(@(x,y,z)sin(x));
F3 = ballfun(@(x,y,z)cos(y.*z));
F = ballfunv(F1,F2,F3);
G = -F;
pass(1) = norm( F + G ) < tol;
end

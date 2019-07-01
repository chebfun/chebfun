function pass = test_power( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Test with function r*cos(lam)*sin(th) to the power 5
f = ballfun(@(r,lam,th)r.*cos(lam).*sin(th), 'spherical');
g = power(f,5);
h = ballfun(@(r,lam,th)r.^5.*cos(lam).^5.*sin(th).^5, 'spherical');

pass(1) = norm( g - h ) < tol;
end

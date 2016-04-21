function pass = test_biharm(pref)
% Check the CHEBFUN2 BIHARM command.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6*pref.cheb2Prefs.chebfun2eps;

% Function to be used:
ff = @(x,y) x.^2.*y.^2;

% Bihamrmonic operator applied to ff:
f = chebfun2(ff); 
fB = biharm(f); 

% Exact solution:
gg = @(x,y) 8;
g = chebfun2(gg);

% Compare:
pass = ( norm(fB - g) < tol );

end
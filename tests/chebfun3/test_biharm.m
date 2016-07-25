function pass = test_biharm(pref)
% Test BIHARM command.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6*pref.cheb3Prefs.chebfun3eps;

% Function to be used:
ff = @(x,y,z) x.^2.*y.^2 + x.^2.*z.^2 + y.^2.*z.^2;

% Bihamrmonic operator applied to ff:
f = chebfun3(ff); 
fB = biharm(f); 

% Exact solution:
gg = @(x,y,z) 24;
g = chebfun3(gg);

% Compare:
pass = norm(fB - g) < tol;

end
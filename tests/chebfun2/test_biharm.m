function pass = test_biharm( pref )
% Check the biharm command in Chebfun2

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e6 * pref.eps; 

% Function to be used
ff = @(x,y) x.^2.*y.^2;

% Bihamrmonic operator applied to ff
gg = @(x,y) 8;
g = chebfun2(gg);

f = chebfun2(ff); 
fB = biharm(f); 

pass = ( norm( fB - g ) < tol);

end